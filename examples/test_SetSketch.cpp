/**
 * test_SetSketch – SetSketch all-to-all with direct pairwise verification.
 *
 * This path intentionally avoids inverted-index candidate pruning to prevent
 * pair leakage. We directly run all-pairs SetSketch Jaccard and report Mash
 * distance.
 *
 * Usage:
 *   exe_test_SetSketch <file_list> <dist_threshold> <threads> [output_file]
 */

#include "Sketch.h"
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
#include <cstring>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <numeric>
#include <sys/stat.h>
#include <unistd.h>

using namespace std;

KSEQ_INIT(gzFile, gzread)

static inline double mash_distance_from_jaccard(double jaccard, int kmerSize) {
    if (jaccard <= 0.0) return 1.0;
    if (jaccard >= 1.0) return 0.0;
    const double p = (2.0 * jaccard) / (1.0 + jaccard);
    return (p > 0.0) ? (-std::log(p) / static_cast<double>(kmerSize)) : 1.0;
}

static inline double min_jaccard_from_mash_distance(double mashDist, int kmerSize) {
    if (mashDist <= 0.0) return 1.0;
    const double p = std::exp(-mashDist * static_cast<double>(kmerSize));
    const double denom = 2.0 - p;
    if (denom <= 0.0) return 1.0;
    double j = p / denom;
    if (j < 0.0) return 0.0;
    if (j > 1.0) return 1.0;
    return j;
}

int main(int argc, char* argv[])
{
    if (argc < 4) {
        cerr << "usage: " << argv[0]
             << " <file_list> <dist_threshold> <threads> [output_file]" << endl;
        return 1;
    }
    const string inputFile = argv[1];
    const double thres     = stod(argv[2]);
    int          nThreads  = stoi(argv[3]);
    if (nThreads < 1) nThreads = 1;
    string outPath = "res.dist.SetSketch";
    if (argc >= 5) outPath = argv[4];

    ifstream fs(inputFile);
    if (!fs) err(errno, "cannot open %s", inputFile.c_str());
    vector<string> fileList;
    { string line; while (getline(fs, line)) if (!line.empty()) fileList.push_back(line); }
    const int N = (int)fileList.size();
    cerr << "===== total files: " << N << "  (SetSketch)" << endl;

    // ── Phase 0: prepare shared sketch parameters + output buffers ───────────
    static const int BITS = 13;
    Sketch::SetSketch proto(BITS, 2.0, 20.0);
    const int m         = proto.getM();
    const double factor = proto.getFactor();
    double bip_buf[64];
    memcpy(bip_buf, proto.getBaseInvPow(), 64 * sizeof(double));
    const double* bip = bip_buf;

    vector<double>   sizes(N, 0.0);
    vector<uint8_t>  flat_cores((size_t)N * m, 0);

    // ── Phase 1a: parallel sketch construction ───────────────────────────────
    double t0 = get_sec();

    #pragma omp parallel for num_threads(nThreads) schedule(dynamic)
    for (int t = 0; t < N; t++) {
        Sketch::SetSketch sk(BITS, 2.0, 20.0);
        gzFile fp = gzopen(fileList[t].c_str(), "r");
        if (!fp) continue;
        kseq_t* ks = kseq_init(fp);
        while (kseq_read(ks) >= 0)
            sk.update(ks->seq.s, static_cast<size_t>(ks->seq.l));
        kseq_destroy(ks);
        gzclose(fp);

        sizes[t] = sk.cardinality();
        memcpy(&flat_cores[(size_t)t * m], sk.getCore().data(), m);
    }

    double t1 = get_sec();
    cerr << "sketch time: " << t1 - t0 << " s" << endl;
    cerr << "flatten + free sketches: skipped (direct build into flat_cores)" << endl;

    const int KMER_SIZE = 32;
    const double minJac = min_jaccard_from_mash_distance(thres, KMER_SIZE);
    const int TILE = 128;

    // Sort by cardinality desc to maximize size-ratio pruning and enable early-break.
    vector<int> order(N);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b) { return sizes[a] > sizes[b]; });
    vector<string> sortedFiles(N);
    vector<double> sortedSizes(N, 0.0);
    vector<uint8_t> sortedCores((size_t)N * m, 0);
    for (int ni = 0; ni < N; ++ni) {
        const int oi = order[ni];
        sortedFiles[ni] = fileList[oi];
        sortedSizes[ni] = sizes[oi];
        memcpy(&sortedCores[(size_t)ni * m], &flat_cores[(size_t)oi * m], m);
    }
    fileList.swap(sortedFiles);
    sizes.swap(sortedSizes);
    flat_cores.swap(sortedCores);

    cerr << "mode: DIRECT all-pairs (no inverted index pruning)\n";
    cerr << "pruning: minJac(from mashDist<" << thres << ") = " << minJac << "\n";

    // Handle output path (directory detection)
    string finalPath = outPath;
    {
        struct stat st;
        if (stat(finalPath.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
            if (finalPath.back() != '/') finalPath += '/';
            finalPath += "res.dist.SetSketch";
        }
    }
    FILE* fout = fopen(finalPath.c_str(), "wb");
    if (!fout) {
        cerr << "ERROR: cannot open output file: " << finalPath << endl;
        return 1;
    }
    setvbuf(fout, nullptr, _IOFBF, 1 << 24);
    cerr << "output: " << finalPath << endl;
    fclose(fout);

    int progress = N / 20;
    if (progress < 1) progress = 1;

    double t5 = get_sec();
    const int nTiles = (N + TILE - 1) / TILE;
    vector<string> partPaths(nThreads);
    const int pid = static_cast<int>(getpid());
    for (int tid = 0; tid < nThreads; ++tid)
        partPaths[tid] = finalPath + ".part." + to_string(pid) + "." + to_string(tid);

    #pragma omp parallel num_threads(nThreads)
    {
        const int tid = omp_get_thread_num();
        FILE* tf = fopen(partPaths[tid].c_str(), "wb");
        if (!tf) {
            #pragma omp critical
            cerr << "ERROR: cannot open temp output file: " << partPaths[tid] << "\n";
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

                    for (int i = iBeg; i < iEnd; ++i) {
                        int jStart = (tileI == tileJ) ? max(i + 1, jBeg) : jBeg;
                        if (jStart >= jEnd) continue;

                        const uint8_t* core_i = &flat_cores[(size_t)i * m];
                        const double si = sizes[i];

                        // sizes sorted desc: in each i-loop, maxJacBySize decreases as j grows.
                        for (int j = jStart; j < jEnd; ++j) {
                            const double sj = sizes[j];
                            const double maxJacBySize = (si > 0.0) ? (sj / si) : 0.0;
                            if (maxJacBySize < minJac) {
                                break;
                            }
                            const uint8_t* c2 = &flat_cores[(size_t)j * m];
                            double jaccard = Sketch::SetSketch::jaccardFromCoresEarlyAbort(
                                core_i, c2, m, bip, factor, si, sj, minJac);
                            if (jaccard < minJac) continue;
                            double dist = mash_distance_from_jaccard(jaccard, KMER_SIZE);

                            if (dist < thres) {
                                char line[1024];
                                int len = snprintf(line, sizeof(line), "%s\t%s\t%.6f\n",
                                                   fileList[i].c_str(), fileList[j].c_str(), dist);
                                buf.append(line, static_cast<size_t>(len));
                            }

                            if (buf.size() > (1 << 22)) {
                                fwrite(buf.data(), 1, buf.size(), tf);
                                buf.clear();
                            }
                        }

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

    fout = fopen(finalPath.c_str(), "wb");
    if (!fout) {
        cerr << "ERROR: cannot open output file for merge: " << finalPath << endl;
        return 1;
    }
    setvbuf(fout, nullptr, _IOFBF, 1 << 24);
    vector<char> copyBuf(1 << 20);
    for (int tid = 0; tid < nThreads; ++tid) {
        FILE* part = fopen(partPaths[tid].c_str(), "rb");
        if (!part) continue;
        while (true) {
            size_t got = fread(copyBuf.data(), 1, copyBuf.size(), part);
            if (got == 0) break;
            fwrite(copyBuf.data(), 1, got, fout);
        }
        fclose(part);
        remove(partPaths[tid].c_str());
    }
    fclose(fout);

    double t6 = get_sec();
    cerr << "dist time: " << t6 - t5 << " s" << endl;
    cerr << "total time: " << t6 - t0 << " s" << endl;

    return 0;
}
