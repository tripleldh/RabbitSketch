/**
 * test_SetSketch – fully optimized benchmark harness.
 *
 * Key optimizations vs. naive version:
 *  1. Pre-allocate sketch array → eliminate omp critical section
 *  2. Flat contiguous core array → cache-friendly band hash + distance
 *  3. Parallel sort for band_entries
 *  4. LSH bucket size cap → prevent candidate explosion
 *  5. Inline SIMD distance on flat array with prefetch
 *  6. schedule(dynamic, chunk) for amortised atomic overhead
 *
 * Usage: exe_test_SetSketch <file_list> <dist_threshold> <threads> [max_bucket=500]
 */
#include "Sketch.h"
#include <sys/time.h>
#include <zlib.h>
#include "kseq.h"
#include <vector>
#include <cmath>
#include <fstream>
#include <err.h>
#include <omp.h>
#include <algorithm>
#include <atomic>
#include <cstring>
#include <immintrin.h>
#include "common.h"
#if defined(__GNUC__) && defined(_OPENMP)
# include <parallel/algorithm>
#endif
using namespace std;

KSEQ_INIT(gzFile, gzread)

// ── Inline union-size on flat core arrays (SIMD gather) ──────────────────────
static inline double flat_union_card(const uint8_t* __restrict__ c1,
                                     const uint8_t* __restrict__ c2,
                                     int m,
                                     const double* __restrict__ bip,
                                     double factor)
{
  double sum = 0.0;
#if defined(__AVX512BW__) && defined(__AVX512F__)
  int i = 0;
  __m512d vsum0 = _mm512_setzero_pd();
  __m512d vsum1 = _mm512_setzero_pd();
  for (; i + 16 <= m; i += 16) {
    __m128i va = _mm_loadu_si128((__m128i*)(c1 + i));
    __m128i vb = _mm_loadu_si128((__m128i*)(c2 + i));
    __m128i vmax = _mm_max_epu8(va, vb);
    __m256i vidx0 = _mm256_cvtepu8_epi32(vmax);
    vsum0 = _mm512_add_pd(vsum0, _mm512_i32gather_pd(vidx0, bip, 8));
    __m128i vmax_hi = _mm_srli_si128(vmax, 8);
    __m256i vidx1 = _mm256_cvtepu8_epi32(vmax_hi);
    vsum1 = _mm512_add_pd(vsum1, _mm512_i32gather_pd(vidx1, bip, 8));
  }
  sum = _mm512_reduce_add_pd(_mm512_add_pd(vsum0, vsum1));
  for (; i < m; i++) {
    uint8_t r = (c1[i] > c2[i]) ? c1[i] : c2[i];
    sum += bip[r];
  }
#else
  for (int i = 0; i < m; i++) {
    uint8_t r = (c1[i] > c2[i]) ? c1[i] : c2[i];
    sum += bip[r];
  }
#endif
  return (sum > 1e-300) ? factor / sum : 0.0;
}

int main(int argc, char* argv[])
{
  if (argc < 4) {
    cerr << "usage: " << argv[0]
         << " <file_list> <dist_threshold> <threads> [max_bucket=500]" << endl;
    return 1;
  }
  string inputFile = argv[1];
  double thres     = stod(argv[2]);
  int numThreads   = stoi(argv[3]);
  if (numThreads < 1) numThreads = 1;
  int MAX_BUCKET   = (argc >= 5) ? stoi(argv[4]) : 500;

  ifstream fs(inputFile);
  if (!fs) err(errno, "cannot open %s", inputFile.c_str());

  vector<string> fileArr;
  { string line; while (getline(fs, line)) fileArr.push_back(line); }
  const int n = (int)fileArr.size();
  cerr << "===== total files: " << n << " (SetSketch, optimized)" << endl;

  // ── Phase 0: pre-allocate sketches (serial, ~0.1 s) ────────────────────────
  static const int BITS = 10;
  vector<Sketch::SetSketch> vsketch;
  vsketch.reserve(n);
  for (int i = 0; i < n; i++)
    vsketch.emplace_back(BITS, 2.0, 20.0);

  // ── Phase 1: parallel sketch construction (NO critical section) ─────────────
  double t1 = get_sec();

  #pragma omp parallel for num_threads(numThreads) schedule(dynamic)
  for (int t = 0; t < n; t++) {
    gzFile fp1 = gzopen(fileArr[t].c_str(), "r");
    if (fp1 == NULL) continue;
    kseq_t* ks1 = kseq_init(fp1);
    while (kseq_read(ks1) >= 0)
      vsketch[t].update(ks1->seq.s);
    kseq_destroy(ks1);
    gzclose(fp1);
  }

  double t2 = get_sec();
  cerr << "sketch time is: " << t2 - t1 << endl;

  // ── Phase 2: extract flat core array + pre-compute cardinalities ────────────
  const int m = vsketch[0].getM();
  const double factor = vsketch[0].getFactor();
  const double* bip   = vsketch[0].getBaseInvPow();

  vector<double> sizes(n);
  vector<uint8_t> flat_cores((size_t)n * m);

  #pragma omp parallel for num_threads(numThreads) schedule(static)
  for (int i = 0; i < n; i++) {
    sizes[i] = vsketch[i].cardinality();
    memcpy(&flat_cores[(size_t)i * m], vsketch[i].getCore().data(), m);
  }
  double t_phase1 = get_sec();
  cerr << "  Phase 1 cardinality + flatten (parallel): " << t_phase1 - t2 << " s" << endl;

  // ── Phase 3: LSH banding on flat array ──────────────────────────────────────
  const int BANDS = 128;
  const int ROWS  = m / BANDS;

  auto band_hash = [](const uint8_t* data, int len) -> uint32_t {
    uint32_t h = 2166136261u;
    for (int i = 0; i < len; i++) { h ^= data[i]; h *= 16777619u; }
    return h;
  };

  vector<pair<uint64_t,int>> band_entries((size_t)n * BANDS);
  #pragma omp parallel for num_threads(numThreads) schedule(static)
  for (int i = 0; i < n; i++) {
    const uint8_t* core = &flat_cores[(size_t)i * m];
    for (int b = 0; b < BANDS; b++) {
      uint32_t h = band_hash(core + b * ROWS, ROWS);
      band_entries[(size_t)i * BANDS + b] =
        { ((uint64_t)(uint32_t)b << 32) | (uint32_t)h, i };
    }
  }
  double t_build = get_sec();
  cerr << "  Phase 2 build band_entries (parallel): " << t_build - t_phase1 << " s" << endl;

  // Parallel sort
#if defined(__GNUC__) && defined(_OPENMP)
  __gnu_parallel::sort(band_entries.begin(), band_entries.end());
#else
  sort(band_entries.begin(), band_entries.end());
#endif
  double t_sort_bands = get_sec();
  cerr << "  Phase 2 sort band_entries [PARALLEL]: " << t_sort_bands - t_build << " s" << endl;

  // Group boundaries
  vector<size_t> group_start;
  group_start.reserve(band_entries.size() / 4);
  group_start.push_back(0);
  for (size_t i = 1; i < band_entries.size(); i++)
    if (band_entries[i].first != band_entries[i-1].first)
      group_start.push_back(i);
  group_start.push_back(band_entries.size());
  const int n_bands = (int)group_start.size() - 1;

  // Generate candidates with bucket cap
  long long skipped_buckets = 0;
  vector<vector<pair<int,int>>> thread_cands((size_t)numThreads);

  #pragma omp parallel num_threads(numThreads)
  {
    int tid = omp_get_thread_num();
    auto& local = thread_cands[tid];
    local.clear();
    #pragma omp for schedule(dynamic) reduction(+:skipped_buckets)
    for (int g = 0; g < n_bands; g++) {
      size_t s = group_start[g], e = group_start[g+1];
      int bsz = (int)(e - s);
      if (bsz > MAX_BUCKET) { skipped_buckets++; continue; }
      for (size_t a = s; a < e; a++)
        for (size_t b = a+1; b < e; b++) {
          int ia = band_entries[a].second, ib = band_entries[b].second;
          local.emplace_back(min(ia, ib), max(ia, ib));
        }
    }
  }
  band_entries.clear();
  band_entries.shrink_to_fit();

  // Merge per-thread candidates
  vector<size_t> prefix((size_t)numThreads + 1);
  prefix[0] = 0;
  for (int t = 0; t < numThreads; t++)
    prefix[t+1] = prefix[t] + thread_cands[t].size();
  size_t total_cands = prefix[numThreads];
  vector<pair<int,int>> candidates(total_cands);
  #pragma omp parallel for num_threads(numThreads)
  for (int t = 0; t < numThreads; t++)
    copy(thread_cands[t].begin(), thread_cands[t].end(),
         candidates.begin() + prefix[t]);
  thread_cands.clear();
  thread_cands.shrink_to_fit();
  double t_scan = get_sec();
  cerr << "  Phase 2 scan->candidates (parallel): " << t_scan - t_sort_bands << " s"
       << "  (skipped " << skipped_buckets << " large buckets)" << endl;

  // Dedup
#if defined(__GNUC__) && defined(_OPENMP)
  __gnu_parallel::sort(candidates.begin(), candidates.end());
#else
  sort(candidates.begin(), candidates.end());
#endif
  candidates.erase(unique(candidates.begin(), candidates.end()), candidates.end());
  double t_lsh = get_sec();
  cerr << "  Phase 2 sort+dedup candidates: " << t_lsh - t_scan << " s" << endl;

  const long long total_pairs = (long long)n * (n - 1) / 2;
  cerr << "LSH candidates: " << candidates.size()
       << " / " << total_pairs << " total pairs"
       << "  (reduction: " << 100.0*(1.0-(double)candidates.size()/total_pairs) << "%)" << endl;
  cerr << "LSH index time (total): " << t_lsh - t2 << " s" << endl;

  // ── Phase 4: thread-local output buffers (zero IO during computation) ─────
  vector<string> thread_bufs(numThreads);

  // ── Phase 5: distance computation (inline SIMD on flat array + prefetch) ───
  const double min_jaccard = 1.0 - thres;
  const int ncand = (int)candidates.size();
  atomic<long long> cnt_size_filtered{0};
  atomic<long long> cnt_exact{0};

  #pragma omp parallel for num_threads(numThreads) schedule(dynamic, 4096)
  for (int c = 0; c < ncand; c++) {
    int i = candidates[c].first, j = candidates[c].second;

    // Size-ratio prefilter
    double si = sizes[i], sj = sizes[j];
    if (si <= 0 || sj <= 0) { cnt_size_filtered++; continue; }
    if (min(si, sj) / max(si, sj) < min_jaccard * 0.93) { cnt_size_filtered++; continue; }

    // Prefetch next pair's core data
    if (c + 1 < ncand) {
      __builtin_prefetch(&flat_cores[(size_t)candidates[c+1].first  * m], 0, 0);
      __builtin_prefetch(&flat_cores[(size_t)candidates[c+1].second * m], 0, 0);
    }

    // Inline distance on flat array
    cnt_exact++;
    const uint8_t* c1 = &flat_cores[(size_t)i * m];
    const uint8_t* c2 = &flat_cores[(size_t)j * m];
    double us = flat_union_card(c1, c2, m, bip, factor);
    if (us <= 0.0) continue;
    double inter = si + sj - us;
    double jaccard = (inter > 0.0) ? inter / us : 0.0;
    double dist = 1.0 - jaccard;

    if (dist < thres) {
      int tid = omp_get_thread_num();
      char line[4096];
      int len = snprintf(line, sizeof(line), "%s\t%s\t%lf\n",
                         fileArr[i].c_str(), fileArr[j].c_str(), dist);
      thread_bufs[tid].append(line, len);
    }
  }

  // Flush all buffers to a single output file
  system("mkdir -p res_dir");
  FILE* fp_out = fopen("res_dir/res.dist.SetSketch", "w");
  for (int t = 0; t < numThreads; t++)
    fwrite(thread_bufs[t].data(), 1, thread_bufs[t].size(), fp_out);
  fclose(fp_out);

  double t3 = get_sec();
  cerr << "dist time is: " << t3 - t_lsh << " s" << endl;
  cerr << "  candidates:    " << (long long)candidates.size() << endl;
  cerr << "  size-filtered: " << cnt_size_filtered.load() << endl;
  cerr << "  exact computed:" << cnt_exact.load() << endl;

  return 0;
}
