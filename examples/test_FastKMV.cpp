/**
 * test_FastKMV – benchmark harness for FastKMV (fast bottom-k KMV sketch).
 *
 * Mirrors test_ProbMinHash.cpp but uses the KMV Jaccard estimator (sorted
 * two-pointer merge) instead of register-wise equality counting.
 *
 * Usage:
 *   exe_test_FastKMV <file_list> <dist_threshold> <threads> [max_bucket=500]
 */

#include "fastkmv.h"
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

// Scalar sorted-set-intersection with budget on distinct elements.
static uint64_t flat_scalar_intersect(const uint64_t* list1, uint64_t size1,
                                      const uint64_t* list2, uint64_t size2,
                                      uint64_t budget,
                                      uint64_t* i_a, uint64_t* i_b)
{
    uint64_t counter = 0;
    *i_a = 0;
    *i_b = 0;
    while (*i_a < size1 && *i_b < size2 && budget > 0) {
        if (list1[*i_a] < list2[*i_b])      { ++(*i_a); }
        else if (list1[*i_a] > list2[*i_b]) { ++(*i_b); }
        else { ++counter; ++(*i_a); ++(*i_b); }
        --budget;
    }
    return counter;
}

// KMV Jaccard on two sorted flat arrays of size k.
// SIMD-accelerated rotational intersection (AVX-512 / AVX2 / scalar).
static inline double flat_kmv_jaccard(const uint64_t* __restrict__ a,
                                      const uint64_t* __restrict__ b,
                                      int k)
{
    const uint64_t sa = static_cast<uint64_t>(
        std::lower_bound(a, a + k, UINT64_MAX) - a);
    const uint64_t sb = static_cast<uint64_t>(
        std::lower_bound(b, b + k, UINT64_MAX) - b);
    const uint64_t K = static_cast<uint64_t>(k);

    if (sa == 0 || sb == 0) return 0.0;

    uint64_t ia = 0, ib = 0;
    uint64_t common = 0;

#if defined(__AVX512F__)
    {
        const uint64_t st_a = (sa / 8) * 8;
        const uint64_t st_b = (sb / 8) * 8;

        if (K > 8 && st_a > 0 && st_b > 0) {
            const uint64_t stop = K - 8;

            __m512i sv0 = _mm512_set_epi64(0,7,6,5,4,3,2,1);
            __m512i sv1 = _mm512_set_epi64(1,0,7,6,5,4,3,2);
            __m512i sv2 = _mm512_set_epi64(2,1,0,7,6,5,4,3);
            __m512i sv3 = _mm512_set_epi64(3,2,1,0,7,6,5,4);
            __m512i sv4 = _mm512_set_epi64(4,3,2,1,0,7,6,5);
            __m512i sv5 = _mm512_set_epi64(5,4,3,2,1,0,7,6);
            __m512i sv6 = _mm512_set_epi64(6,5,4,3,2,1,0,7);

            while (ia < st_a && ib < st_b) {
                __m512i va = _mm512_loadu_si512(&a[ia]);
                __m512i vb = _mm512_loadu_si512(&b[ib]);

                uint64_t a_max = a[ia + 7];
                uint64_t b_max = b[ib + 7];

                ia += (a_max <= b_max) * 8;
                ib += (a_max >= b_max) * 8;

                __mmask8 cmp0 = _mm512_cmpeq_epu64_mask(va, vb);
                __mmask8 cmp1 = _mm512_cmpeq_epu64_mask(va, _mm512_permutexvar_epi64(sv0, vb));
                __mmask8 cmp2 = _mm512_cmpeq_epu64_mask(va, _mm512_permutexvar_epi64(sv1, vb));
                __mmask8 cmp3 = _mm512_cmpeq_epu64_mask(va, _mm512_permutexvar_epi64(sv2, vb));
                __mmask8 cmpA = cmp0 | cmp1 | cmp2 | cmp3;

                __mmask8 cmp4 = _mm512_cmpeq_epu64_mask(va, _mm512_permutexvar_epi64(sv3, vb));
                __mmask8 cmp5 = _mm512_cmpeq_epu64_mask(va, _mm512_permutexvar_epi64(sv4, vb));
                __mmask8 cmp6 = _mm512_cmpeq_epu64_mask(va, _mm512_permutexvar_epi64(sv5, vb));
                __mmask8 cmp7 = _mm512_cmpeq_epu64_mask(va, _mm512_permutexvar_epi64(sv6, vb));
                __mmask8 cmpB = cmp4 | cmp5 | cmp6 | cmp7;

                __mmask8 hits = cmpA | cmpB;
                common += _mm_popcnt_u64(hits);

                if (ia + ib - common >= stop) {
                    common -= _mm_popcnt_u64(hits);
                    ia -= (a_max <= b_max) * 8;
                    ib -= (a_max >= b_max) * 8;
                    break;
                }
            }
        }
    }
#elif defined(__AVX2__)
    {
        const uint64_t st_a = (sa / 4) * 4;
        const uint64_t st_b = (sb / 4) * 4;

        if (K > 8 && st_a > 0 && st_b > 0) {
            const uint64_t stop = K - 8;

            while (ia < st_a && ib < st_b) {
                __m256i va = _mm256_loadu_si256((const __m256i*)&a[ia]);
                __m256i vb = _mm256_loadu_si256((const __m256i*)&b[ib]);

                uint64_t a_max = a[ia + 3];
                uint64_t b_max = b[ib + 3];

                ia += (a_max <= b_max) * 4;
                ib += (a_max >= b_max) * 4;

                __m256i cmp1 = _mm256_cmpeq_epi64(va, vb);
                __m256i rot1 = _mm256_permute4x64_epi64(vb, 0x39);
                __m256i cmp2 = _mm256_cmpeq_epi64(va, rot1);
                __m256i rot2 = _mm256_permute4x64_epi64(vb, 0x4E);
                __m256i cmp3 = _mm256_cmpeq_epi64(va, rot2);
                __m256i rot3 = _mm256_permute4x64_epi64(vb, 0x93);
                __m256i cmp4 = _mm256_cmpeq_epi64(va, rot3);

                __m256i combined = _mm256_or_si256(
                    _mm256_or_si256(cmp1, cmp2),
                    _mm256_or_si256(cmp3, cmp4));
                int mask = _mm256_movemask_pd(_mm256_castsi256_pd(combined));
                common += _mm_popcnt_u64(static_cast<unsigned>(mask));

                if (ia + ib - common >= stop) break;
            }
        }
    }
#endif

    uint64_t remaining = K - (ia + ib - common);
    uint64_t ia_s, ib_s;
    common += flat_scalar_intersect(a + ia, sa - ia, b + ib, sb - ib,
                                    remaining, &ia_s, &ib_s);
    ia += ia_s;
    ib += ib_s;

    uint64_t distinct = ia + ib - common;
    while (distinct < K && ia < sa) { ++distinct; ++ia; }
    while (distinct < K && ib < sb) { ++distinct; ++ib; }
    if (distinct > K) distinct = K;

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
    cerr << "===== total files: " << n << "  (FastKMV)" << endl;

    static const uint32_t K     = 1024;
    static const int      KSIZE = 21;
    static const uint64_t SEED  = 42;
    const int m = (int)K;

    vector<Sketch::FastKMV> vsketches;
    vector<string>          paths;
    vsketches.reserve(n);
    paths.reserve(n);

    double t1 = get_sec();

    #pragma omp parallel for num_threads(numThreads) schedule(dynamic)
    for (int t = 0; t < n; t++) {
        gzFile fp1 = gzopen(fileArr[t].c_str(), "r");
        if (fp1 == NULL) continue;
        kseq_t* ks1 = kseq_init(fp1);

        Sketch::FastKMV sk(K, KSIZE, SEED);
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

    vector<uint64_t> flat_regs((size_t)n_actual * m);

    #pragma omp parallel for num_threads(numThreads) schedule(static)
    for (int i = 0; i < n_actual; i++) {
        const uint64_t* src = vsketches[i].getRegisters();
        memcpy(&flat_regs[(size_t)i * m], src, m * sizeof(uint64_t));
    }
    { vector<Sketch::FastKMV>().swap(vsketches); }

    double t_flat = get_sec();
    cerr << "flatten + free sketches: " << t_flat - t2 << " s" << endl;

    // ── Phase 3: LSH banding ─────────────────────────────────────────────────
    // FNV-1a on sorted value bands.  Less effective for KMV than for
    // register-based sketches, but still generates useful candidates.
    const int BANDS = 128;
    const int ROWS  = m / BANDS;

    auto band_hash = [&](const uint64_t* data, int rows) -> uint32_t {
        uint32_t h = 2166136261u;
        const uint8_t* p = reinterpret_cast<const uint8_t*>(data);
        for (int i = 0; i < rows * (int)sizeof(uint64_t); i++) {
            h ^= p[i];
            h *= 16777619u;
        }
        return h;
    };

    vector<pair<uint64_t,int>> band_entries((size_t)n_actual * BANDS);
    #pragma omp parallel for num_threads(numThreads) schedule(static)
    for (int i = 0; i < n_actual; i++) {
        const uint64_t* regs = &flat_regs[(size_t)i * m];
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
        const uint64_t* ra = &flat_regs[(size_t)i * m];
        const uint64_t* rb = &flat_regs[(size_t)j * m];
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
    FILE* fp_out = fopen("res_dir/res.dist.FastKMV", "w");
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
