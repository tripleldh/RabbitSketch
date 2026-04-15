/**
 * test_BinDash – BinDash all-to-all (library-centric mode, no inverted index).
 *
 * Builds BinDash sketches, then performs direct pairwise calls to library
 * jaccard()/distance logic (SIMD path encapsulated inside src/BinDash.cpp).
 *
 * Usage:
 *   exe_test_BinDash <file_list> <dist_threshold> <threads> [output_file]
 */

#include "BinDash.h"
#include "common.h"
#include "kseq.h"

#include <zlib.h>
#include <sys/time.h>
#include <err.h>
#include <omp.h>

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <optional>

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
    if (N == 0) {
        cerr << "no input files found in: " << inputFile << endl;
        return 1;
    }

    static const uint32_t SKETCH64 = 16;
    static const int      KSIZE    = 21;
    static const uint32_t BBITS    = 16;
    static const uint64_t SEED     = 42;
    const uint32_t NBINS = SKETCH64 * 64;

    const int actualThreads = min(nThreads, N);
    vector<optional<Sketch::BinDash>> sketches(N);

    double t0 = get_sec();

    // ── Phase 1: Sketch + finalize ─────────────────────────────────────────
    #pragma omp parallel num_threads(actualThreads)
    {
        #pragma omp for schedule(dynamic)
        for (int t = 0; t < N; t++) {
            gzFile fp = gzopen(fileList[t].c_str(), "r");
            if (!fp) continue;
            kseq_t* ks = kseq_init(fp);
            sketches[t].emplace(SKETCH64, KSIZE, BBITS, SEED);
            auto& sk = *sketches[t];

            while (kseq_read(ks) >= 0)
                sk.update(ks->seq.s, ks->seq.l);
            sk.finalize();

            kseq_destroy(ks);
            gzclose(fp);
        }
    }

    double t1 = get_sec();
    cerr << "sketch + finalize: " << t1 - t0 << " s" << endl;

    // ── Phase 2: Direct all-to-all distance via library API ───────────────
    const double p_exp  = exp(-static_cast<double>(KSIZE) * maxDist);
    const double minJac = p_exp / (2.0 - p_exp);
    const int minSamebits = max(1,
        static_cast<int>(ceil(NBINS * (minJac * (1.0 - 1.0 / static_cast<double>(1ULL << BBITS))
                                       + 1.0 / static_cast<double>(1ULL << BBITS)))));
    cerr << "pruning: minJac=" << minJac << "  minSamebits=" << minSamebits
         << "/" << NBINS << "  (mashD<" << maxDist << ", k=" << KSIZE << ")" << endl;

    FILE* fout = fopen(outPath.c_str(), "w");
    if (!fout) err(errno, "cannot open %s", outPath.c_str());
    setvbuf(fout, nullptr, _IOFBF, 1 << 24);

    double t2 = get_sec();
    #pragma omp parallel num_threads(actualThreads)
    {
        string buf;
        buf.reserve(1 << 16); // 64KB per thread, lower peak memory than 1MB/thread

        #pragma omp for schedule(dynamic, 8)
        for (int i = 0; i < N; ++i) {
            if (!sketches[i]) continue;
            const auto& si = *sketches[i];
            for (int j = i + 1; j < N; ++j) {
                if (!sketches[j]) continue;
                const auto& sj = *sketches[j];
                const double jac = si.jaccard(sj);
                if (jac < minJac) continue;

                const double dist = (jac >= 1.0) ? 0.0
                    : -std::log(2.0 * jac / (1.0 + jac)) / static_cast<double>(KSIZE);
                if (dist < maxDist) {
                    char line[1024];
                    int n = snprintf(line, sizeof(line), "%s\t%s\t%.6f\n",
                                     fileList[i].c_str(), fileList[j].c_str(), dist);
                    buf.append(line, static_cast<size_t>(n));
                }
            }
            if (buf.size() > (1 << 16)) {
                #pragma omp critical
                { fwrite(buf.data(), 1, buf.size(), fout); }
                buf.clear();
            }
        }
        if (!buf.empty()) {
            #pragma omp critical
            { fwrite(buf.data(), 1, buf.size(), fout); }
        }
    }
    fclose(fout);
    double t3 = get_sec();

    cerr << "dist time: " << t3 - t2 << " s" << endl;
    cerr << "total time: " << t3 - t0 << " s" << endl;
    cerr << "output: " << outPath << endl;

    return 0;
}
