/**
 * test_BinDash – BinDash all-to-all with flat sketch storage.
 *
 * Memory layout (vs old vector<optional<BinDash>>):
 *   Old: N separate heap allocs (usigs_) + N × ~81-byte struct overhead
 *   New: one contiguous flat_usigs array + N × 9 bytes (size + valid flag)
 *   Savings: N × ~81 bytes (80 MB at N=1M, 800 MB at N=10M)
 *
 * Additional optimizations:
 *   - Sort by cardinality → break inner loop when size ratio < minJac
 *   - Tile-level cardinality pruning (skip whole tile blocks)
 *   - minSamebits pre-check before Jaccard formula (skips rejected pairs faster)
 *   - Per-thread temp files → single fwrite bottleneck eliminated
 *   - Prefetch next sketch into cache
 *
 * Usage:
 *   exe_test_BinDash <file_list> <dist_threshold> <threads> [output_file]
 */

#include "BinDash.h"
#include "common.h"
#include "kseq.h"

#include <zlib.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <err.h>
#include <omp.h>
#include <unistd.h>
#include <signal.h>

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>

// ── Global state for temp-file cleanup on signals ──────────────────────────
// SIGTERM / SIGINT: remove in-flight part files before dying.
static std::vector<std::string>* g_partPaths = nullptr;
static void sig_cleanup(int) {
    if (g_partPaths)
        for (const auto& p : *g_partPaths) remove(p.c_str());
    _exit(1);
}

using namespace std;

KSEQ_INIT(gzFile, gzread)

int main(int argc, char* argv[])
{
    if (argc < 4) {
        cerr << "usage: " << argv[0]
             << " <file_list> <dist_threshold> <threads> [output_file]" << endl;
        return 1;
    }
    const string inputFile = argv[1];
    const double maxDist   = stod(argv[2]);
    int          nThreads  = stoi(argv[3]);
    if (nThreads < 1) nThreads = 1;
    const string outPath = (argc >= 5) ? argv[4] : "res.dist.BinDash";

    ifstream fs(inputFile);
    if (!fs) err(errno, "cannot open %s", inputFile.c_str());
    vector<string> fileList;
    { string line; while (getline(fs, line)) if (!line.empty()) fileList.push_back(line); }
    const int N = (int)fileList.size();
    cerr << "===== total files: " << N << "  (BinDash)" << endl;
    if (N == 0) { cerr << "no input files found\n"; return 1; }

    // ── Sketch parameters ────────────────────────────────────────────────────
    static const uint32_t SKETCH64 = 16;
    static const int      KSIZE    = 21;
    static const uint32_t BBITS    = 16;
    static const uint64_t SEED     = 42;
    static const uint32_t NBINS    = SKETCH64 * 64;    // 1024 bins
    static const uint32_t NWORDS   = SKETCH64 * BBITS; // 256 uint64_t per sketch

    const int actualThreads = min(nThreads, N);

    // ── Flat sketch storage ─────────────────────────────────────────────────
    // One contiguous block instead of N separate heap-allocated usigs_ vectors.
    // Eliminates N allocator calls and N × ~81-byte BinDash struct overhead.
    vector<uint64_t> flat_usigs((size_t)N * NWORDS, 0);
    vector<double>   sizes(N, 0.0);  // cardinality estimate for size-ratio pruning
    vector<uint8_t>  valid(N, 0);

    double t0 = get_sec();

    // ── Phase 1: Sketch + extract into flat array ───────────────────────────
    #pragma omp parallel for num_threads(actualThreads) schedule(dynamic)
    for (int t = 0; t < N; t++) {
        gzFile fp = gzopen(fileList[t].c_str(), "r");
        if (!fp) continue;
        kseq_t* ks = kseq_init(fp);

        Sketch::BinDash sk(SKETCH64, KSIZE, BBITS, SEED);
        while (kseq_read(ks) >= 0)
            sk.update(ks->seq.s, ks->seq.l);
        sk.finalize();  // densify + pack + free signs_ (8 KB temp)

        kseq_destroy(ks);
        gzclose(fp);

        memcpy(&flat_usigs[(size_t)t * NWORDS], sk.getSignatures(),
               NWORDS * sizeof(uint64_t));

        // Cardinality from occupied-bin count (same formula as BinDash::cardinalityEstimate)
        const double ne = static_cast<double>(sk.getRawNonempty());
        sizes[t] = (ne <= 0.0) ? 0.0
                 : (ne >= NBINS) ? 1e18
                 : -static_cast<double>(NBINS) * std::log1p(-ne / NBINS);
        valid[t] = 1;
        // sk goes out of scope here; no BinDash object persists
    }

    double t1 = get_sec();
    cerr << "sketch + finalize: " << t1 - t0 << " s" << endl;

    // ── Sort by cardinality descending ──────────────────────────────────────
    // After sorting: sizes[i] >= sizes[j] for i < j.
    // Inner j loop can break when sizes[j] / sizes[i] < minJac (impossible to meet
    // Jaccard threshold given containment bound min(|A|,|B|)/max(|A|,|B|) < minJac).
    vector<int> order(N);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b) { return sizes[a] > sizes[b]; });

    // In-place permutation via cycle following.
    // Uses only one 2048-byte sketch buffer instead of a full N × NWORDS copy.
    // For N=208K saves ~428 MB of peak RSS vs the copy approach.
    {
        vector<uint64_t> tmp(NWORDS);
        vector<bool> done(N, false);
        for (int i = 0; i < N; ++i) {
            if (done[i] || order[i] == i) { done[i] = true; continue; }
            // Save element i
            memcpy(tmp.data(), &flat_usigs[(size_t)i * NWORDS], NWORDS * 8);
            double   tmp_sz = sizes[i];
            string   tmp_f  = std::move(fileList[i]);
            uint8_t  tmp_v  = valid[i];
            int j = i;
            int k = order[j];
            while (k != i) {
                memcpy(&flat_usigs[(size_t)j * NWORDS],
                       &flat_usigs[(size_t)k * NWORDS], NWORDS * 8);
                sizes[j]    = sizes[k];
                fileList[j] = std::move(fileList[k]);
                valid[j]    = valid[k];
                done[j]     = true;
                j = k;
                k = order[j];
            }
            memcpy(&flat_usigs[(size_t)j * NWORDS], tmp.data(), NWORDS * 8);
            sizes[j]    = tmp_sz;
            fileList[j] = std::move(tmp_f);
            valid[j]    = tmp_v;
            done[j]     = true;
        }
    }
    { vector<int>().swap(order); }  // free 4×N bytes, no longer needed

    // ── Resolve output path (handle directory) ─────────────────────────────
    string finalPath = outPath;
    {
        struct stat st;
        if (stat(finalPath.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
            if (finalPath.back() != '/') finalPath += '/';
            finalPath += "res.dist.BinDash";
        }
    }

    // ── Pruning thresholds ──────────────────────────────────────────────────
    const double p_exp   = std::exp(-static_cast<double>(KSIZE) * maxDist);
    const double minJac  = p_exp / (2.0 - p_exp);
    const double p_rand  = 1.0 / static_cast<double>(1ULL << BBITS);  // 1/65536
    // samebits >= minSamebits ⟺ jac >= minJac (exact with ceiling)
    const uint64_t minSamebits = static_cast<uint64_t>(std::max(1,
        static_cast<int>(std::ceil(NBINS * (minJac * (1.0 - p_rand) + p_rand)))));
    cerr << "pruning: minJac=" << minJac
         << "  minSamebits=" << minSamebits << "/" << NBINS
         << "  (mashD<" << maxDist << ", k=" << KSIZE << ")" << endl;

    // ── Per-thread temp output files ────────────────────────────────────────
    const int procId = static_cast<int>(getpid());
    vector<string> partPaths(actualThreads);
    for (int tid = 0; tid < actualThreads; ++tid)
        partPaths[tid] = finalPath + ".part." + to_string(procId) + "." + to_string(tid);

    // Register signal handlers so temp files are cleaned up on SIGTERM/SIGINT
    g_partPaths = &partPaths;
    signal(SIGTERM, sig_cleanup);
    signal(SIGINT,  sig_cleanup);

    double t2 = get_sec();

    // Tile size: 64 sketches × 2048 bytes = 128 KB → fits comfortably in L2
    static const int TILE = 64;
    const int nTiles = (N + TILE - 1) / TILE;
    int progress = N / 20;
    if (progress < 1) progress = 1;

    // ── Phase 2: Tiled distance computation ────────────────────────────────
    #pragma omp parallel num_threads(actualThreads)
    {
        const int tid = omp_get_thread_num();
        FILE* tf = fopen(partPaths[tid].c_str(), "wb");
        if (!tf) {
            #pragma omp critical
            cerr << "ERROR: cannot open temp file: " << partPaths[tid] << "\n";
        } else {
            setvbuf(tf, nullptr, _IOFBF, 1 << 22);
            string buf;
            buf.reserve(1 << 22);

            #pragma omp for schedule(dynamic, 1)
            for (int tileI = 0; tileI < nTiles; ++tileI) {
                const int iBeg = tileI * TILE;
                const int iEnd = min(N, iBeg + TILE);

                for (int tileJ = tileI; tileJ < nTiles; ++tileJ) {
                    const int jBeg = tileJ * TILE;
                    const int jEnd = min(N, jBeg + TILE);

                    // Tile-level prune: largest in jBeg vs smallest in iEnd-1
                    if (tileJ > tileI && sizes[iEnd - 1] > 0.0 &&
                        sizes[jBeg] / sizes[iEnd - 1] < minJac)
                        break;

                    for (int i = iBeg; i < iEnd; ++i) {
                        if (!valid[i]) continue;
                        const uint64_t* ui = &flat_usigs[(size_t)i * NWORDS];
                        const double    si = sizes[i];

                        const int jStart = (tileI == tileJ) ? max(i + 1, jBeg) : jBeg;

                        for (int j = jStart; j < jEnd; ++j) {
                            if (!valid[j]) continue;

                            // Per-pair size-ratio abort
                            if (si > 0.0 && sizes[j] / si < minJac) break;

                            if (j + 4 < jEnd)
                                __builtin_prefetch(&flat_usigs[(size_t)(j + 4) * NWORDS], 0, 1);

                            const uint64_t* uj = &flat_usigs[(size_t)j * NWORDS];

                            // Count matching b-bit values (SIMD XNOR + popcount)
                            const uint64_t same = Sketch::BinDash::countSameBits(
                                ui, uj, SKETCH64, BBITS);
                            if (same < minSamebits) continue;

                            // Jaccard from samebits (same >= minSamebits ⟹ jac >= minJac)
                            const double p_match = static_cast<double>(same) / NBINS;
                            const double jac     = (p_match - p_rand) / (1.0 - p_rand);
                            const double dist    = (jac >= 1.0) ? 0.0
                                : -std::log(2.0 * jac / (1.0 + jac))
                                  / static_cast<double>(KSIZE);

                            if (dist < maxDist) {
                                char line[1024];
                                int len = snprintf(line, sizeof(line), "%s\t%s\t%.6f\n",
                                    fileList[i].c_str(), fileList[j].c_str(), dist);
                                buf.append(line, static_cast<size_t>(len));
                            }
                        }
                    }

                    if (buf.size() > (1 << 22)) {
                        fwrite(buf.data(), 1, buf.size(), tf);
                        buf.clear();
                    }
                }

                const int iProgress = min(iEnd, N);
                if (iProgress % progress == 0) {
                    #pragma omp critical
                    cerr << "  dist " << iProgress << " / " << N << "\n";
                }
            }
            if (!buf.empty()) fwrite(buf.data(), 1, buf.size(), tf);
            fclose(tf);
        }
    }

    // ── Merge temp files ─────────────────────────────────────────────────────
    // Always clean up part files before returning, even on error.
    auto cleanup_parts = [&]() {
        for (int tid = 0; tid < actualThreads; ++tid) remove(partPaths[tid].c_str());
    };

    FILE* fout = fopen(finalPath.c_str(), "wb");
    if (!fout) {
        cerr << "ERROR: cannot open output: " << finalPath << "\n";
        cleanup_parts();
        return 1;
    }
    setvbuf(fout, nullptr, _IOFBF, 1 << 24);
    vector<char> copyBuf(1 << 20);
    for (int tid = 0; tid < actualThreads; ++tid) {
        FILE* part = fopen(partPaths[tid].c_str(), "rb");
        if (!part) continue;
        while (true) {
            size_t got = fread(copyBuf.data(), 1, copyBuf.size(), part);
            if (!got) break;
            fwrite(copyBuf.data(), 1, got, fout);
        }
        fclose(part);
        remove(partPaths[tid].c_str());
    }
    fclose(fout);
    g_partPaths = nullptr;  // disarm signal handler

    double t3 = get_sec();
    cerr << "dist time: " << t3 - t2 << " s" << endl;
    cerr << "total time: " << t3 - t0 << " s" << endl;
    cerr << "output: " << finalPath << endl;

    return 0;
}
