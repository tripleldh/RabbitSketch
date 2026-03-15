/**
 * test_ProbMinHash – benchmark harness for ProbMinHash4 sketch.
 *
 * Mirrors the structure of test_HLL.cpp and test_SetSketch.cpp:
 *   Phase 0 : pre-allocate sketch objects
 *   Phase 1 : parallel sketch construction (rolling canonical k-mers)
 *   Phase 2 : LSH banding on flat register arrays → candidate pairs
 *   Phase 3 : verify candidates with exact distance()
 *
 * Probability Jaccard distance = 1 - J_p.
 * The Jaccard estimator is: fraction of registers that carry identical
 * minimum-hash values (standard ProbMinHash bottom-k estimator).
 *
 * Usage:
 *   exe_test_ProbMinHash <file_list> <dist_threshold> <threads> [max_bucket=500]
 *
 * <file_list>      : one genome file path per line (gzipped FASTA/FASTQ)
 * <dist_threshold> : output pairs with distance < threshold (e.g. 0.1)
 * <threads>        : number of OpenMP threads
 * [max_bucket]     : LSH bucket cap to avoid candidate explosion (default 500)
 */

#include "probmh.h"
#include "common.h"
#include "kseq.h"

#include <zlib.h>
#include <sys/time.h>
#include <err.h>
#include <omp.h>

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <atomic>
#include <cstring>

#if defined(__GNUC__) && defined(_OPENMP)
#  include <parallel/algorithm>
#endif

using namespace std;

KSEQ_INIT(gzFile, gzread)

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
    cerr << "===== total files: " << n << "  (ProbMinHash4)" << endl;

    // ── Sketch parameters ────────────────────────────────────────────────────
    static const uint32_t M     = 1024;  // number of registers
    static const int      KSIZE = 21;    // k-mer length
    static const uint64_t SEED  = 42;

    // ── Phase 0: pre-allocate sketch objects (serial) ────────────────────────
    vector<Sketch::ProbMinHash4> vsketches;
    vsketches.reserve(n);
    for (int i = 0; i < n; i++)
        vsketches.emplace_back(M, KSIZE, SEED);

    // ── Phase 1: parallel sketch construction ────────────────────────────────
    double t1 = get_sec();

    #pragma omp parallel for num_threads(numThreads) schedule(dynamic)
    for (int t = 0; t < n; t++) {
        gzFile  fp1 = gzopen(fileArr[t].c_str(), "r");
        if (fp1 == NULL) continue;
        kseq_t* ks1 = kseq_init(fp1);
        while (kseq_read(ks1) >= 0)
            vsketches[t].update(ks1->seq.s);
        kseq_destroy(ks1);
        gzclose(fp1);
    }

    double t2 = get_sec();
    cerr << "sketch time: " << t2 - t1 << " s" << endl;

    // ── Phase 2: extract flat register arrays for cache-friendly LSH ─────────
    // Each sketch stores M doubles; pack them row-major into a flat array.
    const int m = (int)M;
    vector<double> flat_regs((size_t)n * m);

    #pragma omp parallel for num_threads(numThreads) schedule(static)
    for (int i = 0; i < n; i++) {
        const double* src = vsketches[i].getRegisters().data();
        memcpy(&flat_regs[(size_t)i * m], src, m * sizeof(double));
    }
    double t_flat = get_sec();
    cerr << "  flatten registers (parallel): " << t_flat - t2 << " s" << endl;

    // ── Phase 3: LSH banding ─────────────────────────────────────────────────
    // Hash each band of ROWS consecutive register values into a bucket key.
    // Two sketches that share at least one band hash are candidate pairs.
    const int BANDS = 128;
    const int ROWS  = m / BANDS;  // 8 registers per band with M=1024, BANDS=128

    // Band hash: FNV-1a over the raw double bytes of one band
    auto band_hash = [&](const double* data, int rows) -> uint32_t {
        uint32_t h = 2166136261u;
        const uint8_t* p = reinterpret_cast<const uint8_t*>(data);
        for (int i = 0; i < rows * (int)sizeof(double); i++) {
            h ^= p[i];
            h *= 16777619u;
        }
        return h;
    };

    vector<pair<uint64_t,int>> band_entries((size_t)n * BANDS);
    #pragma omp parallel for num_threads(numThreads) schedule(static)
    for (int i = 0; i < n; i++) {
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
    cerr << "  sort band_entries: " << t_sort - t_build << " s" << endl;

    // Group boundaries
    vector<size_t> group_start;
    group_start.reserve(band_entries.size() / 4);
    group_start.push_back(0);
    for (size_t i = 1; i < band_entries.size(); i++)
        if (band_entries[i].first != band_entries[i-1].first)
            group_start.push_back(i);
    group_start.push_back(band_entries.size());
    const int n_groups = (int)group_start.size() - 1;

    // Generate candidate pairs with bucket cap
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

    // Merge per-thread candidates into contiguous array
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

    // Dedup
#if defined(__GNUC__) && defined(_OPENMP)
    __gnu_parallel::sort(candidates.begin(), candidates.end());
#else
    sort(candidates.begin(), candidates.end());
#endif
    candidates.erase(unique(candidates.begin(), candidates.end()), candidates.end());
    double t_lsh = get_sec();
    cerr << "  sort+dedup candidates: " << t_lsh - t_scan << " s" << endl;

    const long long total_pairs = (long long)n * (n - 1) / 2;
    cerr << "LSH candidates: " << candidates.size()
         << " / " << total_pairs << " total pairs"
         << "  (reduction: "
         << 100.0 * (1.0 - (double)candidates.size() / (double)total_pairs)
         << "%)" << endl;
    cerr << "LSH index time (total): " << t_lsh - t2 << " s" << endl;

    // ── Phase 4: verify candidates with exact distance ───────────────────────
    // Thread-local string buffers → zero disk IO during computation.
    vector<string> thread_bufs(numThreads);
    const int    ncand       = (int)candidates.size();
    atomic<long long> cnt_exact{0};

    #pragma omp parallel for num_threads(numThreads) schedule(dynamic, 4096)
    for (int c = 0; c < ncand; c++) {
        int i = candidates[c].first, j = candidates[c].second;
        int tid = omp_get_thread_num();

        cnt_exact++;
        double dist = vsketches[i].distance(vsketches[j]);
        if (dist < thres) {
            char line[4096];
            int len = snprintf(line, sizeof(line), "%s\t%s\t%.6f\n",
                               fileArr[i].c_str(), fileArr[j].c_str(), dist);
            thread_bufs[tid].append(line, len);
        }
    }

    // Flush all buffers to a single output file
    system("mkdir -p res_dir");
    FILE* fp_out = fopen("res_dir/res.dist.ProbMinHash", "w");
    for (int t = 0; t < numThreads; t++)
        fwrite(thread_bufs[t].data(), 1, thread_bufs[t].size(), fp_out);
    fclose(fp_out);

    double t3 = get_sec();
    cerr << "dist time: " << t3 - t_lsh << " s" << endl;
    cerr << "  candidates:      " << (long long)candidates.size() << endl;
    cerr << "  output pairs:    " << cnt_exact.load() << endl;
    cerr << "total time: " << t3 - t1 << " s" << endl;

    return 0;
}
