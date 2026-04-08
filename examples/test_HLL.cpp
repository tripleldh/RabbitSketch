#include "Sketch.h"
#include <sys/time.h>
#include <zlib.h>
#include "kseq.h"
#include <vector>
#include <cmath>
#include <fstream>
#include <err.h>
#include <sys/stat.h>
#include <omp.h>
#include <atomic>
#include <algorithm>
#include <cstring>
#include "common.h"
#include "robin_hood.h"
using namespace std;

KSEQ_INIT(gzFile, gzread)


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
  cerr << "===== total files: " << n << " (HyperLogLog, optimized)" << endl;

  // ── Phase 0: pre-allocate sketches ──────────────────────────────────────────
  static const int BITS = 13;
  vector<Sketch::HyperLogLog> vhlog;
  vhlog.reserve(n);
  for (int i = 0; i < n; i++)
    vhlog.emplace_back(BITS);

  // ── Phase 1: parallel sketch construction (no critical section) ─────────────
  double t1 = get_sec();

  #pragma omp parallel for num_threads(numThreads) schedule(dynamic)
  for (int t = 0; t < n; t++) {
    gzFile fp1 = gzopen(fileArr[t].c_str(), "r");
    if (fp1 == NULL) continue;
    kseq_t* ks1 = kseq_init(fp1);
    while (kseq_read(ks1) >= 0)
      vhlog[t].update(ks1->seq.s);
    kseq_destroy(ks1);
    gzclose(fp1);
  }

  double t2 = get_sec();
  cerr << "sketch time is: " << t2 - t1 << endl;

  // ── Phase 2: extract flat core array + cardinalities ────────────────────────
  const int m = (n > 0) ? (int)vhlog[0].getCore().size() : 0;
  vector<double> sizes(n);
  vector<uint8_t> flat_cores((size_t)n * m);

  #pragma omp parallel for num_threads(numThreads) schedule(static)
  for (int i = 0; i < n; i++) {
    sizes[i] = vhlog[i].cardinality();
    memcpy(&flat_cores[(size_t)i * m], vhlog[i].getCore().data(), m);
  }
  double t_phase1 = get_sec();
  cerr << "  Phase 1 cardinality + flatten (parallel): " << t_phase1 - t2 << " s" << endl;

  // ── Phase 3–5: Band-streaming LSH + inline verification ───────────────────
  const int BANDS = 128;
  const int ROWS  = (m > 0) ? m / BANDS : 8;
  const double min_jaccard = 1.0 - thres;

  auto band_hash = [](const uint8_t* data, int len) -> uint32_t {
    uint32_t h = 2166136261u;
    for (int i = 0; i < len; i++) { h ^= data[i]; h *= 16777619u; }
    return h;
  };

  vector<string> thread_bufs(numThreads);

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
  long long skipped_buckets = 0;

  double t_lsh_start = get_sec();

  for (int b = 0; b < BANDS; b++) {
    robin_hood::unordered_map<uint32_t, vector<int>> bkt;
    bkt.reserve((size_t)n);
    for (int i = 0; i < n; i++) {
      uint32_t h = band_hash(flat_cores.data() + (size_t)i * m + b * ROWS, ROWS);
      bkt[h].push_back(i);
    }

    vector<pair<int,int>> new_pairs;
    for (auto& [key, ids] : bkt) {
      int sz = (int)ids.size();
      if (sz < 2) continue;
      if (sz > MAX_BUCKET) { skipped_buckets++; continue; }
      for (int a = 0; a < sz; a++)
        for (int c = a + 1; c < sz; c++) {
          int ii = min(ids[a], ids[c]), jj = max(ids[a], ids[c]);
          double si = sizes[ii], sj = sizes[jj];
          if (si <= 0 || sj <= 0) continue;
          if (min(si, sj) / max(si, sj) < min_jaccard * 0.93) {
            cnt_size_filtered++;
            continue;
          }
          new_pairs.emplace_back(ii, jj);
        }
    }

    const int np = (int)new_pairs.size();
    #pragma omp parallel for num_threads(numThreads) schedule(dynamic, 256)
    for (int c = 0; c < np; c++) {
      int i = new_pairs[c].first, j = new_pairs[c].second;
      if (!try_insert(i, j)) continue;

      cnt_exact++;
      double dist = vhlog[i].distance(vhlog[j]);
      if (dist < thres) {
        int tid = omp_get_thread_num();
        char line[4096];
        int len = snprintf(line, sizeof(line), "%s\t%s\t%lf\n",
                           fileArr[i].c_str(), fileArr[j].c_str(), dist);
        thread_bufs[tid].append(line, len);
      }
    }
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

  // Flush output
  system("mkdir -p res_dir");
  FILE* fp_out = fopen("res_dir/res.dist.HLL", "w");
  for (int t = 0; t < numThreads; t++)
    fwrite(thread_bufs[t].data(), 1, thread_bufs[t].size(), fp_out);
  fclose(fp_out);

  cerr << "total time: " << get_sec() - t1 << " s" << endl;
  return 0;
}




















