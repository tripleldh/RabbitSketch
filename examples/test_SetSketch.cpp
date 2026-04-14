/**
 * test_SetSketch – SetSketch all-to-all via block-of-3 inverted index.
 *
 * SetSketch registers are 8-bit (only 256 values), so individual registers
 * lack entropy for effective inverted-index filtering.  We group every 3
 * adjacent registers into a block and hash (block_idx, v1, v2, v3) into a
 * uint32_t key.  Random collision probability drops to ~(1/256)^3 ≈ 6e-8,
 * while similar genomes (J≈0.2) still share ~25 matching blocks on average.
 *
 * Candidates passing the block-match threshold are verified with exact
 * SetSketch Jaccard — zero accuracy loss.
 *
 * Usage:
 *   exe_test_SetSketch <file_list> <dist_threshold> <threads> [output_file]
 */

#include "Sketch.h"
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
#include <cstring>
#include <cmath>
#include <climits>
#include <cstdint>
#include <sys/stat.h>

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
    const double thres     = stod(argv[2]);
    int          nThreads  = stoi(argv[3]);
    if (nThreads < 1) nThreads = 1;
    const string outPath = (argc >= 5) ? argv[4] : "res.dist.SetSketch";

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

    vector<double>  sizes(N, 0.0);
    vector<uint8_t> flat_cores((size_t)N * m, 0);

    // ── Phase 1a: parallel sketch construction ───────────────────────────────
    double t0 = get_sec();

    #pragma omp parallel for num_threads(nThreads) schedule(dynamic)
    for (int t = 0; t < N; t++) {
        Sketch::SetSketch sk(BITS, 2.0, 20.0);
        gzFile fp = gzopen(fileList[t].c_str(), "r");
        if (!fp) continue;
        kseq_t* ks = kseq_init(fp);
        while (kseq_read(ks) >= 0)
            sk.update(ks->seq.s);
        kseq_destroy(ks);
        gzclose(fp);

        sizes[t] = sk.cardinality();
        memcpy(&flat_cores[(size_t)t * m], sk.getCore().data(), m);
    }

    double t1 = get_sec();
    cerr << "sketch time: " << t1 - t0 << " s" << endl;
    double t2 = t1;
    cerr << "flatten + free sketches: skipped (direct build into flat_cores)" << endl;

    // ── Phase 1c: build local inverted index with block-of-3 keys ────────────
    const int NUM_BLOCKS = Sketch::SetSketch::numBlocks(m);

    const int actualThreads = min(nThreads, N);
    vector<phmap::flat_hash_map<uint32_t, vector<uint32_t>>> threadIdx(actualThreads);

    #pragma omp parallel num_threads(actualThreads)
    {
        int tid = omp_get_thread_num();
        auto& localIdx = threadIdx[tid];

        #pragma omp for schedule(static)
        for (int t = 0; t < N; t++) {
            const uint8_t* core = &flat_cores[(size_t)t * m];
            for (int b = 0; b < NUM_BLOCKS; b++) {
                uint32_t key = Sketch::SetSketch::blockHash(
                    static_cast<uint32_t>(b),
                    core[b * 3], core[b * 3 + 1], core[b * 3 + 2]);
                localIdx[key].push_back(static_cast<uint32_t>(t));
            }
        }
    }

    double t3 = get_sec();
    cerr << "local inverted index: " << t3 - t2 << " s" << endl;

    // ── Phase 2: Build CSR inverted index ────────────────────────────────────
    auto csrIdx = Sketch::buildCSRIndex<uint32_t>(threadIdx, nThreads);
    double t4 = get_sec();
    cerr << "build CSR index: " << t4 - t3 << " s" << endl;

    // ── Phase 3: distance with exact verification ────────────────────────────
    const int minMatchBlocks = Sketch::SetSketch::minMatchBlocksForDist(thres, NUM_BLOCKS);

    cerr << "pruning: minMatchBlocks=" << minMatchBlocks << "/" << NUM_BLOCKS
         << "  (minJac=" << (1.0 - thres) << ", maxDist=" << thres << ")" << endl;

    // Handle output path (directory detection)
    string finalPath = outPath;
    {
        struct stat st;
        if (stat(finalPath.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
            if (finalPath.back() != '/') finalPath += '/';
            finalPath += "res.dist.SetSketch";
        }
    }
    FILE* fout = fopen(finalPath.c_str(), "w");
    if (!fout) {
        cerr << "ERROR: cannot open output file: " << finalPath << endl;
        return 1;
    }
    setvbuf(fout, nullptr, _IOFBF, 1 << 24);
    cerr << "output: " << finalPath << endl;

    const uint32_t* csrPtr = csrIdx.csrPosts.data();
    int progress = N / 20;
    if (progress < 1) progress = 1;

    double t5 = get_sec();

    #pragma omp parallel num_threads(nThreads)
    {
        vector<uint16_t> isect(N, 0);
        vector<int> stamp(N, 0);
        int ep = 0;
        vector<int> cand;
        cand.reserve(4096);
        string buf;
        buf.reserve(1 << 24);

        #pragma omp for schedule(dynamic, 64)
        for (int i = 0; i < N; i++) {
            cand.clear();
            ++ep;
            if (__builtin_expect(ep == INT_MAX, 0)) {
                memset(stamp.data(), 0, N * sizeof(int));
                ep = 1;
            }

            // Reconstruct block keys from flat_cores (no skKeys storage needed)
            const uint8_t* core_i = &flat_cores[(size_t)i * m];
            for (int b = 0; b < NUM_BLOCKS; b++) {
                uint32_t key = Sketch::SetSketch::blockHash(
                    static_cast<uint32_t>(b),
                    core_i[b * 3], core_i[b * 3 + 1], core_i[b * 3 + 2]);

                auto it = csrIdx.postIdx.find(key);
                if (__builtin_expect(it == csrIdx.postIdx.end(), 0)) continue;

                const uint32_t* pl   = csrPtr + it->second.off;
                const uint32_t  plSz = it->second.cnt;
                for (uint32_t pi = 0; pi < plSz; pi++) {
                    int j = static_cast<int>(pl[pi]);
                    if (j <= i) continue;
                    if (__builtin_expect(stamp[j] != ep, 1)) {
                        stamp[j] = ep;
                        isect[j] = 1;
                        cand.push_back(j);
                    } else {
                        isect[j]++;
                    }
                }
            }

            // Verify candidates with exact SetSketch Jaccard
            const double si = sizes[i];
            for (int j : cand) {
                if (isect[j] < minMatchBlocks) continue;

                const uint8_t* c2 = &flat_cores[(size_t)j * m];
                double jaccard = Sketch::SetSketch::jaccardFromCores(
                    core_i, c2, m, bip, factor, si, sizes[j]);
                double dist = 1.0 - jaccard;

                if (dist < thres) {
                    char line[1024];
                    int len = snprintf(line, sizeof(line), "%s\t%s\t%.6f\n",
                                       fileList[i].c_str(), fileList[j].c_str(), dist);
                    buf.append(line, static_cast<size_t>(len));
                }
            }

            if (buf.size() > (1 << 24)) {
                #pragma omp critical
                { fwrite(buf.data(), 1, buf.size(), fout); }
                buf.clear();
            }
            if (i % progress == 0)
                cerr << "  dist " << i << " / " << N << "\n";
        }
        if (!buf.empty()) {
            #pragma omp critical
            { fwrite(buf.data(), 1, buf.size(), fout); }
        }
    }
    fclose(fout);

    double t6 = get_sec();
    cerr << "dist time: " << t6 - t5 << " s" << endl;
    cerr << "total time: " << t6 - t0 << " s" << endl;

    return 0;
}
