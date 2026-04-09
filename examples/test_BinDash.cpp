/**
 * test_BinDash – BinDash all-to-all via inverted index (zero false negatives).
 *
 * Builds an inverted index from each bin's truncated bbits-bit value to the
 * genomes that share it.  Key = (bin_idx << 16) | 16bit_value.
 * Distance is computed by counting matching bins (samebits) per pair via
 * posting-list traversal, then applying the same b-bit correction as
 * BinDash::jaccard().
 *
 * Usage:
 *   exe_test_BinDash <file_list> <dist_threshold> <threads> [output_file]
 */

#include "BinDash.h"
#include "InvertedIndex.h"
#include "common.h"
#include "kseq.h"
#include "phmap.h"

#include <zlib.h>
#include <sys/time.h>
#include <err.h>
#include <omp.h>

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

KSEQ_INIT(gzFile, gzread)

int main(int argc, char* argv[])
{
    if (argc < 4) {
        cerr << "usage: " << argv[0]
             << " <file_list> <dist_threshold> <threads> [output_file]" << endl;
        return 1;
    }
    const string inputFile  = argv[1];
    const double maxDist    = stod(argv[2]);
    int          nThreads   = stoi(argv[3]);
    if (nThreads < 1) nThreads = 1;
    const string outPath = (argc >= 5) ? argv[4] : "res.dist.BinDash";

    ifstream fs(inputFile);
    if (!fs) err(errno, "cannot open %s", inputFile.c_str());
    vector<string> fileList;
    { string line; while (getline(fs, line)) if (!line.empty()) fileList.push_back(line); }
    const int N = (int)fileList.size();
    cerr << "===== total files: " << N << "  (BinDash)" << endl;

    static const uint32_t SKETCH64 = 16;
    static const int      KSIZE    = 21;
    static const uint32_t BBITS    = 16;
    static const uint64_t SEED     = 42;
    const uint32_t NBINS = SKETCH64 * 64;

    // Per-genome signature: extracted bin values for inverted-index lookup
    vector<vector<uint32_t>> skKeys(N);
    vector<int> sketchSizes(N, static_cast<int>(NBINS));

    const int actualThreads = min(nThreads, N);
    vector<phmap::flat_hash_map<uint32_t, vector<uint32_t>>> threadIdx(actualThreads);

    double t0 = get_sec();

    // ── Phase 1: Sketch + finalize + local inverted index ─────────────────
    #pragma omp parallel num_threads(nThreads)
    {
        int tid = omp_get_thread_num();
        auto& localIdx = threadIdx[tid];

        #pragma omp for schedule(dynamic)
        for (int t = 0; t < N; t++) {
            gzFile fp = gzopen(fileList[t].c_str(), "r");
            if (!fp) continue;
            kseq_t* ks = kseq_init(fp);

            Sketch::BinDash sk(SKETCH64, KSIZE, BBITS, SEED);
            while (kseq_read(ks) >= 0)
                sk.update(ks->seq.s, ks->seq.l);
            sk.finalize();

            kseq_destroy(ks);
            gzclose(fp);

            skKeys[t].resize(NBINS);
            for (uint32_t bin = 0; bin < NBINS; bin++) {
                uint32_t key = (bin << 16) | sk.getBinValue(bin);
                skKeys[t][bin] = key;
                localIdx[key].push_back(static_cast<uint32_t>(t));
            }
        }
    }

    double t1 = get_sec();
    cerr << "sketch + local index: " << t1 - t0 << " s" << endl;

    // ── Phase 2: Build CSR inverted index ────────────────────────────────
    auto csrIdx = Sketch::buildCSRIndex<uint32_t>(threadIdx, nThreads);
    double t2 = get_sec();
    cerr << "build CSR index: " << t2 - t1 << " s" << endl;

    // ── Phase 3: Distance via inverted index ─────────────────────────────
    const double p_random   = 1.0 / static_cast<double>(1ULL << BBITS);
    const double maxnbits_d = static_cast<double>(NBINS);
    // Derive minimum Jaccard from Mash distance threshold, then minimum samebits
    const double p_exp  = exp(-static_cast<double>(KSIZE) * maxDist);
    const double minJac = p_exp / (2.0 - p_exp);
    // J = (samebits/NBINS - p_random) / (1 - p_random)  ≥  minJac
    // => samebits ≥ NBINS * (minJac*(1-p_random) + p_random)
    const int minSamebits = max(1,
        static_cast<int>(ceil(NBINS * (minJac * (1.0 - p_random) + p_random))));
    cerr << "pruning: minJac=" << minJac << "  minSamebits=" << minSamebits
         << "/" << NBINS << "  (mashD<" << maxDist << ", k=" << KSIZE << ")" << endl;

    auto bbitJaccard = [p_random, maxnbits_d](int samebits, int /*s0*/, int /*s1*/) -> double {
        double p_match = static_cast<double>(samebits) / maxnbits_d;
        return (p_match - p_random) / (1.0 - p_random);
    };

    auto minCommonFn = [minSamebits](int /*s0*/) -> int {
        return minSamebits;
    };

    double t3 = get_sec();
    Sketch::computeDistances<uint32_t>(
        csrIdx, skKeys, sketchSizes, fileList,
        N, KSIZE, maxDist, bbitJaccard, minCommonFn, outPath, nThreads);
    double t4 = get_sec();

    cerr << "dist time: " << t4 - t3 << " s" << endl;
    cerr << "total time: " << t4 - t0 << " s" << endl;
    cerr << "output: " << outPath << endl;

    return 0;
}
