/**
 * test_SetSketch – fully optimized benchmark harness.
 *
 * Key optimizations:
 *  1. Pre-allocate sketch array → eliminate omp critical section
 *  2. Flat contiguous core array → cache-friendly band hash + distance
 *  3. Band-streaming LSH: process one band at a time with hash-map grouping
 *     → eliminates O(n·B) band_entries array + expensive global sort
 *  4. robin_hood hash set for O(1) candidate dedup
 *     → eliminates O(C·log C) candidate sort + dedup
 *  5. Inline SIMD distance on flat array with prefetch
 *  6. LSH bucket size cap → prevent candidate explosion
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
#include "robin_hood.h"
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

  // ── Phase 3–5: Band-streaming LSH + inline verification ─────────────────────
  //
  // Instead of materializing n×BANDS entries, sorting, deduplicating, and then
  // verifying, we process ONE BAND AT A TIME:
  //   1. Build hash-map buckets for band b   → O(n)
  //   2. Enumerate candidate pairs per bucket → O(pairs_in_band)
  //   3. Dedup via robin_hood hash set        → O(1) per pair
  //   4. Verify new pairs inline (parallel)   → O(new_pairs × m)
  //
  // Memory savings:
  //   - Eliminates band_entries  (was n×128×12 bytes)
  //   - Eliminates candidates    (was up to billions of entries + sort)
  //   - Dedup set grows only to unique candidate count

  const int BANDS = 128;
  const int ROWS  = m / BANDS;
  const double min_jaccard = 1.0 - thres;

  auto band_hash = [](const uint8_t* data, int len) -> uint32_t {
    uint32_t h = 2166136261u;
    for (int i = 0; i < len; i++) { h ^= data[i]; h *= 16777619u; }
    return h;
  };

  vector<string> thread_bufs(numThreads);

  // Partitioned dedup: split pair space into NUM_PARTS buckets, each with its
  // own hash set + spinlock, to allow high concurrent insert throughput.
  constexpr int NUM_PARTS = 256;
  struct alignas(64) DedupPart {
    robin_hood::unordered_set<uint64_t> set;
    omp_lock_t lock;
  };
  vector<DedupPart> dedup(NUM_PARTS);
  for (auto& d : dedup) omp_init_lock(&d.lock);

  auto pair_key = [](int i, int j) -> uint64_t {
    return ((uint64_t)(unsigned)i << 32) | (unsigned)j;
  };
  auto try_insert = [&](int i, int j) -> bool {
    uint64_t k = pair_key(i, j);
    int part = (int)((k * 0x9E3779B97F4A7C15ULL) >> 56) & (NUM_PARTS - 1);
    omp_set_lock(&dedup[part].lock);
    bool inserted = dedup[part].set.insert(k).second;
    omp_unset_lock(&dedup[part].lock);
    return inserted;
  };

  atomic<long long> cnt_size_filtered{0};
  atomic<long long> cnt_exact{0};
  long long total_unique_cands = 0;
  long long skipped_buckets = 0;

  double t_lsh_start = get_sec();

  for (int b = 0; b < BANDS; b++) {
    // ── 3a. Build bucket map for this band ──────────────────────────────────
    robin_hood::unordered_map<uint32_t, vector<int>> bkt;
    bkt.reserve((size_t)n);
    for (int i = 0; i < n; i++) {
      uint32_t h = band_hash(flat_cores.data() + (size_t)i * m + b * ROWS, ROWS);
      bkt[h].push_back(i);
    }

    // ── 3b. Collect new candidate pairs (dedup against previous bands) ──────
    vector<pair<int,int>> new_pairs;
    for (auto& [key, ids] : bkt) {
      int sz = (int)ids.size();
      if (sz < 2) continue;
      if (sz > MAX_BUCKET) { skipped_buckets++; continue; }
      for (int a = 0; a < sz; a++)
        for (int c = a + 1; c < sz; c++) {
          int ii = min(ids[a], ids[c]), jj = max(ids[a], ids[c]);
          // Size-ratio prefilter (cheap, avoids hash set lookup)
          double si = sizes[ii], sj = sizes[jj];
          if (si <= 0 || sj <= 0) continue;
          if (min(si, sj) / max(si, sj) < min_jaccard * 0.93) {
            cnt_size_filtered++;
            continue;
          }
          new_pairs.emplace_back(ii, jj);
        }
    }

    // ── 3c. Dedup + inline verification (parallel) ──────────────────────────
    const int np = (int)new_pairs.size();

    #pragma omp parallel for num_threads(numThreads) schedule(dynamic, 256)
    for (int c = 0; c < np; c++) {
      int i = new_pairs[c].first, j = new_pairs[c].second;

      if (!try_insert(i, j)) continue;  // already processed in a previous band

      // Prefetch next pair's core data
      if (c + 1 < np) {
        __builtin_prefetch(&flat_cores[(size_t)new_pairs[c+1].first  * m], 0, 0);
        __builtin_prefetch(&flat_cores[(size_t)new_pairs[c+1].second * m], 0, 0);
      }

      cnt_exact++;
      const uint8_t* c1 = &flat_cores[(size_t)i * m];
      const uint8_t* c2 = &flat_cores[(size_t)j * m];
      double us = flat_union_card(c1, c2, m, bip, factor);
      if (us <= 0.0) continue;
      double si = sizes[i], sj = sizes[j];
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
    total_unique_cands += np;
  }

  for (auto& d : dedup) omp_destroy_lock(&d.lock);
  long long dedup_total = 0;
  for (auto& d : dedup) dedup_total += (long long)d.set.size();

  double t_lsh = get_sec();
  const long long total_pairs = (long long)n * (n - 1) / 2;
  cerr << "LSH band-streaming done (" << BANDS << " bands):" << endl;
  cerr << "  unique candidates:  " << dedup_total << " / " << total_pairs
       << " total pairs (reduction: "
       << 100.0*(1.0-(double)dedup_total/total_pairs) << "%)" << endl;
  cerr << "  skipped buckets:    " << skipped_buckets << endl;
  cerr << "  size-filtered:      " << cnt_size_filtered.load() << endl;
  cerr << "  exact computed:     " << cnt_exact.load() << endl;
  cerr << "  LSH + verify time:  " << t_lsh - t_phase1 << " s" << endl;

  // Flush all buffers to a single output file
  system("mkdir -p res_dir");
  FILE* fp_out = fopen("res_dir/res.dist.SetSketch", "w");
  for (int t = 0; t < numThreads; t++)
    fwrite(thread_bufs[t].data(), 1, thread_bufs[t].size(), fp_out);
  fclose(fp_out);

  double t3 = get_sec();
  cerr << "total time: " << t3 - t1 << " s" << endl;

  return 0;
}
