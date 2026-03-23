/**
 * test_ProbKMV – benchmark harness for ProbKMV (Route B).
 *
 * Mirrors test_ProbMinHash.cpp but uses the KMV Jaccard estimator (sorted
 * two-pointer merge) instead of register-wise equality counting.
 *
 * Usage:
 *   exe_test_ProbKMV <file_list> <dist_threshold> <threads> [max_bucket=500]
 */

#include "probmh.h"
#include "common.h"
#include "kseq.h"

#include <zlib.h>
#include <sys/time.h>
#include <err.h>
#include <omp.h>
#include <immintrin.h>

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <atomic>
#include <cstring>
#include <limits>

#if defined(__GNUC__) && defined(_OPENMP)
#  include <parallel/algorithm>
#endif

using namespace std;

KSEQ_INIT(gzFile, gzread)

// KMV Jaccard on two sorted flat arrays of size k.
// Two-pointer merge: walk the k smallest distinct values from the union,
// count how many appear in both.
static inline double flat_kmv_jaccard(const double* __restrict__ a,
                                      const double* __restrict__ b,
                                      int k)
{
    const double inf = numeric_limits<double>::infinity();
    int ia = 0, ib = 0;
    int distinct = 0, common = 0;

    while (distinct < k && ia < k && ib < k) {
        if (a[ia] >= inf && b[ib] >= inf) break;
        if (a[ia] == b[ib]) {
            ++common; ++distinct; ++ia; ++ib;
        } else if (a[ia] < b[ib]) {
            ++distinct; ++ia;
        } else {
            ++distinct; ++ib;
        }
    }
    while (distinct < k && ia < k && a[ia] < inf) { ++distinct; ++ia; }
    while (distinct < k && ib < k && b[ib] < inf) { ++distinct; ++ib; }

    return (distinct > 0) ? (double)common / (double)distinct : 0.0;
}

int main(int argc, char* argv[])
{
    if (argc < 4) {
        cerr << "usage: " << argv[0]
             << " <file_list> <dist_threshold> <threads> [max_bucket=500]"
             << endl;
        return 1;
    }
    const string inputFile  = argv[1];
    const double thres      = stod(argv[2]);
    int          numThreads = stoi(argv[3]);
    if (numThreads < 1) numThreads = 1;
    const int MAX_BUCKET = (argc >= 5) ? stoi(argv[4]) : 500;

    ifstream fs(inputFile);
    if (!fs) err(errno, "cannot open %s", inputFile.c_str());

    vector<string> fileArr;
    { string line; while (getline(fs, line)) if (!line.empty()) fileArr.push_back(line); }
    const int n = (int)fileArr.size();
    cerr << "===== total files: " << n << "  (ProbKMV – Route B)" << endl;

    static const uint32_t K     = 1024;
    static const int      KSIZE = 21;
    static const uint64_t SEED  = 42;
    const int m = (int)K;

    vector<Sketch::ProbKMV> vsketches;
    vector<string>          paths;
    vsketches.reserve(n);
    paths.reserve(n);

    double t1 = get_sec();

    #pragma omp parallel for num_threads(numThreads) schedule(dynamic)
    for (int t = 0; t < n; t++) {
        gzFile fp1 = gzopen(fileArr[t].c_str(), "r");
        if (fp1 == NULL) continue;
        kseq_t* ks1 = kseq_init(fp1);

        Sketch::ProbKMV sk(K, KSIZE, SEED);
        while (kseq_read(ks1) >= 0)
            sk.update(ks1->seq.s, ks1->seq.l);

        #pragma omp critical
        {
            vsketches.push_back(std::move(sk));
            paths.push_back(fileArr[t]);
        }
        kseq_destroy(ks1);
        gzclose(fp1);
    }

    double t2 = get_sec();
    cerr << "sketch time: " << t2 - t1 << " s" << endl;

    const int n_actual = (int)vsketches.size();
    if (n_actual == 0) { cerr << "no sketches built" << endl; return 1; }

    {
        uint32_t total_fill = 0, min_fill = K;
        for (int i = 0; i < n_actual; i++) {
            uint32_t s = vsketches[i].size();
            total_fill += s;
            if (s < min_fill) min_fill = s;
        }
        cerr << "fill level: avg=" << (double)total_fill / n_actual
             << "  min=" << min_fill << " / " << K << endl;
    }

    vector<double> flat_regs((size_t)n_actual * m);

    #pragma omp parallel for num_threads(numThreads) schedule(static)
    for (int i = 0; i < n_actual; i++) {
        const double* src = vsketches[i].getRegisters();
        memcpy(&flat_regs[(size_t)i * m], src, m * sizeof(double));
    }
    { vector<Sketch::ProbKMV>().swap(vsketches); }

    double t_flat = get_sec();
    cerr << "flatten + free sketches: " << t_flat - t2 << " s" << endl;

    // ── Phase 3: LSH banding ─────────────────────────────────────────────────
    // FNV-1a on sorted value bands.  Less effective for KMV than for
    // register-based sketches, but still generates useful candidates.
    const int BANDS = 128;
    const int ROWS  = m / BANDS;

    auto band_hash = [&](const double* data, int rows) -> uint32_t {
        uint32_t h = 2166136261u;
        const uint8_t* p = reinterpret_cast<const uint8_t*>(data);
        for (int i = 0; i < rows * (int)sizeof(double); i++) {
            h ^= p[i];
            h *= 16777619u;
        }
        return h;
    };

    vector<pair<uint64_t,int>> band_entries((size_t)n_actual * BANDS);
    #pragma omp parallel for num_threads(numThreads) schedule(static)
    for (int i = 0; i < n_actual; i++) {
        const double* regs = &flat_regs[(size_t)i * m];
        for (int b = 0; b < BANDS; b++) {
            uint32_t h = band_hash(regs + b * ROWS, ROWS);
            band_entries[(size_t)i * BANDS + b] =
                { ((uint64_t)(uint32_t)b << 32) | (uint32_t)h, i };
        }
    }
    double t_build = get_sec();
    cerr << "  build band_entries (parallel): " << t_build - t_flat << " s" << endl;

#if defined(__GNUC__) && defined(_OPENMP)
    __gnu_parallel::sort(band_entries.begin(), band_entries.end());
#else
    sort(band_entries.begin(), band_entries.end());
#endif
    double t_sort = get_sec();
    cerr << "  sort band_entries [PARALLEL]: " << t_sort - t_build << " s" << endl;

    vector<size_t> group_start;
    group_start.reserve(band_entries.size() / 4);
    group_start.push_back(0);
    for (size_t i = 1; i < band_entries.size(); i++)
        if (band_entries[i].first != band_entries[i-1].first)
            group_start.push_back(i);
    group_start.push_back(band_entries.size());
    const int n_groups = (int)group_start.size() - 1;

    long long skipped_buckets = 0;
    vector<vector<pair<int,int>>> thread_cands((size_t)numThreads);

    #pragma omp parallel num_threads(numThreads)
    {
        int tid = omp_get_thread_num();
        auto& local = thread_cands[tid];
        local.clear();
        #pragma omp for schedule(dynamic) reduction(+:skipped_buckets)
        for (int g = 0; g < n_groups; g++) {
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
    cerr << "  scan->candidates (parallel): " << t_scan - t_sort << " s"
         << "  (skipped " << skipped_buckets << " large buckets)" << endl;

#if defined(__GNUC__) && defined(_OPENMP)
    __gnu_parallel::sort(candidates.begin(), candidates.end());
#else
    sort(candidates.begin(), candidates.end());
#endif
    candidates.erase(unique(candidates.begin(), candidates.end()), candidates.end());
    double t_lsh = get_sec();
    cerr << "  sort+dedup candidates: " << t_lsh - t_scan << " s" << endl;

    const long long total_pairs = (long long)n_actual * (n_actual - 1) / 2;
    cerr << "LSH candidates: " << candidates.size()
         << " / " << total_pairs << " total pairs"
         << "  (reduction: "
         << 100.0 * (1.0 - (double)candidates.size() / (double)total_pairs)
         << "%)" << endl;
    cerr << "LSH index time (total): " << t_lsh - t2 << " s" << endl;

    // ── Phase 4: verify candidates with KMV Jaccard ──────────────────────────
    vector<string> thread_bufs(numThreads);
    const int    ncand = (int)candidates.size();
    atomic<long long> cnt_exact{0};

    #pragma omp parallel for num_threads(numThreads) schedule(dynamic, 4096)
    for (int c = 0; c < ncand; c++) {
        int i = candidates[c].first, j = candidates[c].second;
        int tid = omp_get_thread_num();

        if (c + 1 < ncand) {
            __builtin_prefetch(&flat_regs[(size_t)candidates[c+1].first  * m], 0, 1);
            __builtin_prefetch(&flat_regs[(size_t)candidates[c+1].second * m], 0, 1);
        }

        cnt_exact++;
        const double* ra = &flat_regs[(size_t)i * m];
        const double* rb = &flat_regs[(size_t)j * m];
        double jaccard = flat_kmv_jaccard(ra, rb, m);
        double dist    = 1.0 - jaccard;

        if (dist < thres) {
            char line[4096];
            int len = snprintf(line, sizeof(line), "%s\t%s\t%.6f\n",
                               paths[i].c_str(), paths[j].c_str(), dist);
            thread_bufs[tid].append(line, len);
        }
    }

    system("mkdir -p res_dir");
    FILE* fp_out = fopen("res_dir/res.dist.ProbKMV", "w");
    for (int t = 0; t < numThreads; t++)
        fwrite(thread_bufs[t].data(), 1, thread_bufs[t].size(), fp_out);
    fclose(fp_out);

    double t3 = get_sec();
    cerr << "dist time: " << t3 - t_lsh << " s" << endl;
    cerr << "  candidates:   " << (long long)candidates.size() << endl;
    cerr << "  exact computed: " << cnt_exact.load() << endl;
    cerr << "total time: " << t3 - t1 << " s" << endl;

    return 0;
}
