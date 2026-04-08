#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <climits>
#include <omp.h>
#include <sys/stat.h>
#ifdef __x86_64__
#include <immintrin.h>
#endif
#include "kseq.h"
#include "zlib.h"
#include "Sketch.h"
#include "phmap.h"
#include "shuffle.h"
#include "common.h"

KSEQ_INIT(gzFile, gzread);

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <file.list> <threshold> <threads> [output_file]\n";
        return 1;
    }
    const std::string inputFile = argv[1];
    const double      maxDist   = std::stod(argv[2]);
    const int         nThreads  = std::stoi(argv[3]);
    const std::string outPath   = (argc >= 5) ? argv[4]
                                              : "result.sketch.dist";

    // ── KSSD parameters ──────────────────────────────────────────────────
    Sketch::kssd_parameter_t P;
    const int  kmerSize = P.kmer_size;
    const bool use64    = (P.half_k - P.drlevel) > 8;
    const uint64_t tupmask = P.tupmask, domask = P.domask;
    const uint64_t umask0 = P.undomask0, umask1 = P.undomask1;
    const int revmov = P.rev_add_move;
    const int domS   = P.half_outctx_len * 2;
    const int drS    = P.drlevel * 4;
    const int undoS  = kmerSize * 2 - P.half_outctx_len * 4;
    const int dimStart = P.dim_start;
    const auto& smap   = P.shuffled_map;

    std::cerr << "half_k=" << P.half_k << " half_subk=" << P.half_subk
              << " drlevel=" << P.drlevel
              << " kmer=" << kmerSize << " use64=" << use64 << "\n";

    // Scalar encode table: ACGT/acgt → 0-3, else 0xFF
    uint8_t BM_U8[256];
    std::memset(BM_U8, 0xFF, sizeof(BM_U8));
    BM_U8['A'] = BM_U8['a'] = 0;
    BM_U8['C'] = BM_U8['c'] = 1;
    BM_U8['G'] = BM_U8['g'] = 2;
    BM_U8['T'] = BM_U8['t'] = 3;

    // ── Read file list ───────────────────────────────────────────────────
    std::vector<std::string> fileList;
    {
        std::ifstream ifs(inputFile);
        if (!ifs) { std::cerr << "ERROR: cannot open " << inputFile << "\n"; return 1; }
        std::string line;
        while (std::getline(ifs, line)) fileList.push_back(line);
    }
    const int N = static_cast<int>(fileList.size());
    std::cerr << "===== total files: " << N << "\n";

    // ── Phase 1: Sketch ─────────────────────────────────────────────────
    // Each genome stores a sorted hash list (uint32 or uint64).
    struct GSketch {
        std::vector<uint32_t> h32;
        std::vector<uint64_t> h64;
    };
    std::vector<GSketch> sketches(N);

    double t0 = get_sec();
    int progress = N / 20; if (progress < 1) progress = 1;

    // Per-thread local inverted indices, merged after the parallel section
    const int actualThreads = std::min(nThreads, N);
    std::vector<phmap::flat_hash_map<uint64_t, std::vector<uint32_t>>>
        threadIdx(actualThreads);

#pragma omp parallel num_threads(nThreads)
    {
        const int tid = omp_get_thread_num();
        phmap::flat_hash_set<uint32_t> hs32;
        phmap::flat_hash_set<uint64_t> hs64;
        hs32.reserve(1 << 16);  // pre-allocate ~64K slots
        hs64.reserve(1 << 16);
        auto& localIdx = threadIdx[tid];

        // SIMD constants: hoist to thread scope (avoid per-sequence reinit)
#if defined(__AVX512BW__)
        const __m512i vlut512 = _mm512_broadcast_i32x4(_mm_setr_epi8(
            -1, 0, -1, 1, 3, -1, -1, 2,
            -1,-1, -1,-1,-1, -1, -1,-1));
        const __m512i vmask512 = _mm512_set1_epi8(0x0F);
#elif defined(__AVX2__)
        const __m256i vlut256 = _mm256_broadcastsi128_si256(_mm_setr_epi8(
            -1, 0, -1, 1, 3, -1, -1, 2,
            -1,-1, -1,-1,-1, -1, -1,-1));
        const __m256i vmask256 = _mm256_set1_epi8(0x0F);
#endif

        // Thread-local I/O buffers for plain FASTA fast path
        std::vector<char> fileBuf;
        std::vector<char> seqBuf;

        // Shared k-mer chunk processor (avoids code duplication)
        auto processSeq = [&](const char* s, int len) {
            if (len < kmerSize) return;
            alignas(64) uint8_t enc[4096];
            uint64_t fwd = 0, rev = 0;
            int run = 0;
            for (int cs = 0; cs < len; cs += 4096) {
                const int clen = std::min(4096, len - cs);
                const char* cp = s + cs;
                int p = 0;
#if defined(__AVX512BW__)
                for (; p + 64 <= clen; p += 64) {
                    __m512i v = _mm512_loadu_si512(cp + p);
                    v = _mm512_shuffle_epi8(vlut512,
                            _mm512_and_si512(v, vmask512));
                    _mm512_storeu_si512(enc + p, v);
                }
#elif defined(__AVX2__)
                for (; p + 32 <= clen; p += 32) {
                    __m256i v = _mm256_loadu_si256(
                            (const __m256i*)(cp + p));
                    v = _mm256_shuffle_epi8(vlut256,
                            _mm256_and_si256(v, vmask256));
                    _mm256_storeu_si256((__m256i*)(enc + p), v);
                }
#endif
                for (; p < clen; ++p)
                    enc[p] = BM_U8[(unsigned char)cp[p]];
                for (int i = 0; i < clen; i++) {
                    const uint8_t b = enc[i];
                    if (__builtin_expect(b <= 3, 1)) {
                        fwd = ((fwd << 2) | b) & tupmask;
                        rev = (rev >> 2)
                            + (((uint64_t)(b ^ 3)) << revmov);
                        if (__builtin_expect(++run >= kmerSize, 1)) {
                            const uint64_t u =
                                (fwd < rev) ? fwd : rev;
                            const uint32_t d = static_cast<uint32_t>(
                                (u & domask) >> domS);
                            auto sit = smap.find(d);
                            if (__builtin_expect(
                                    sit == smap.end(), 1))
                                continue;
                            const uint64_t h =
                                (((u & umask0)
                                  | ((u & umask1) << undoS))
                                 >> drS)
                                | static_cast<uint64_t>(
                                      sit->second - dimStart);
                            if (use64) hs64.emplace(h);
                            else hs32.emplace(
                                     static_cast<uint32_t>(h));
                        }
                    } else {
                        run = 0; fwd = 0; rev = 0;
                    }
                }
            }
        };

#pragma omp for schedule(dynamic)
        for (int t = 0; t < N; t++) {
            hs32.clear(); hs64.clear();

            // Read entire file with fread, detect gzip by magic bytes
            struct stat st;
            if (stat(fileList[t].c_str(), &st) != 0) continue;
            const size_t fsize = static_cast<size_t>(st.st_size);
            if (fsize < 2) continue;
            fileBuf.resize(fsize);
            FILE* ff = fopen(fileList[t].c_str(), "rb");
            if (!ff) continue;
            const size_t rd = fread(fileBuf.data(), 1, fsize, ff);
            fclose(ff);

            const bool isGz =
                (static_cast<unsigned char>(fileBuf[0]) == 0x1f
              && static_cast<unsigned char>(fileBuf[1]) == 0x8b);

            if (isGz) {
                // Gzipped: fall back to kseq streaming reader
                gzFile fp = gzopen(fileList[t].c_str(), "r");
                if (!fp) continue;
                gzbuffer(fp, 1 << 18);
                kseq_t* ks = kseq_init(fp);
                while (kseq_read(ks) >= 0)
                    processSeq(ks->seq.s, static_cast<int>(ks->seq.l));
                kseq_destroy(ks);
                gzclose(fp);
            } else {
                // Plain FASTA: parse directly from memory buffer
                const char* data = fileBuf.data();
                size_t pos = 0;
                while (pos < rd) {
                    if (data[pos] != '>') { pos++; continue; }
                    while (pos < rd && data[pos] != '\n') pos++;
                    if (pos < rd) pos++;
                    seqBuf.clear();
                    while (pos < rd && data[pos] != '>') {
                        size_t ls = pos;
                        while (pos < rd && data[pos] != '\n'
                               && data[pos] != '\r') pos++;
                        seqBuf.insert(seqBuf.end(),
                                      data + ls, data + pos);
                        while (pos < rd && (data[pos] == '\n'
                               || data[pos] == '\r')) pos++;
                    }
                    processSeq(seqBuf.data(),
                               static_cast<int>(seqBuf.size()));
                }
            }

            // Hash set → sketch + inverted index in one pass (no sort needed)
            if (use64) {
                sketches[t].h64.reserve(hs64.size());
                for (uint64_t h : hs64) {
                    sketches[t].h64.push_back(h);
                    localIdx[h].push_back(static_cast<uint32_t>(t));
                }
            } else {
                sketches[t].h32.reserve(hs32.size());
                for (uint32_t h : hs32) {
                    sketches[t].h32.push_back(h);
                    localIdx[static_cast<uint64_t>(h)].push_back(
                        static_cast<uint32_t>(t));
                }
            }
            if (t % progress == 0)
                std::cerr << "  sketch " << t << " / " << N << "\n";
        }
    }
    double t1 = get_sec();
    std::cerr << "sketch + local index: " << t1 - t0 << " s\n";

    // ── Parallel merge via sharded inverted index ──────────────────────
    //   64 shards, each with its own omp_lock → 64 threads can merge
    //   concurrently with minimal contention.
    double t2 = get_sec();
    const int NUM_SHARDS = 64;
    const uint64_t SHARD_MASK = NUM_SHARDS - 1;
    std::vector<phmap::flat_hash_map<uint64_t, std::vector<uint32_t>>>
        invShards(NUM_SHARDS);
    {
        omp_lock_t locks[64];
        for (int s = 0; s < NUM_SHARDS; s++) omp_init_lock(&locks[s]);

#pragma omp parallel num_threads(nThreads)
        {
            int tid = omp_get_thread_num();
            for (auto& [key, vec] : threadIdx[tid]) {
                int shard = static_cast<int>(key & SHARD_MASK);
                omp_set_lock(&locks[shard]);
                auto [it, ins] = invShards[shard].try_emplace(key, std::move(vec));
                if (!ins)
                    it->second.insert(it->second.end(), vec.begin(), vec.end());
                omp_unset_lock(&locks[shard]);
            }
            phmap::flat_hash_map<uint64_t, std::vector<uint32_t>>()
                .swap(threadIdx[tid]);
        }
        for (int s = 0; s < NUM_SHARDS; s++) omp_destroy_lock(&locks[s]);
        threadIdx.clear();
    }
    size_t totalUnique = 0;
    for (auto& sh : invShards) totalUnique += sh.size();
    std::cerr << "merge index: " << totalUnique << " unique hashes, "
              << get_sec() - t2 << " s\n";

    // Parallel singleton removal (one shard per thread, no contention)
    size_t totalBefore = totalUnique, totalAfter = 0;
#pragma omp parallel for num_threads(nThreads) reduction(+:totalAfter)
    for (int s = 0; s < NUM_SHARDS; s++) {
        for (auto it = invShards[s].begin(); it != invShards[s].end(); ) {
            if (it->second.size() <= 1) it = invShards[s].erase(it);
            else { ++it; totalAfter++; }
        }
    }
    std::cerr << "singleton removal: " << totalBefore << " → " << totalAfter
              << " (" << (totalBefore - totalAfter) << " removed), "
              << get_sec() - t2 << " s\n";

    // ── Flatten posting lists into contiguous CSR array ──────────────
    //   RabbitKSSD's key insight: all posting lists packed in one flat
    //   array eliminates 5M scattered heap allocations → eliminates
    //   TLB misses + DRAM random access during dist traversal.
    double tFlat = get_sec();
    size_t totalPostings = 0;
    for (int s = 0; s < NUM_SHARDS; s++)
        for (auto& [k, v] : invShards[s])
            totalPostings += v.size();
    std::cerr << "total non-singleton postings: " << totalPostings
              << " (" << totalPostings * 4.0 / (1ULL << 30) << " GB)\n";

    // Direct array indexing for non-use64: hash value as array index
    // like RabbitKSSD's sketchSizeArr[hash] / offset[hash]
    const uint32_t HASH_SPACE = use64 ? 0
        : (1u << std::min(4 * (P.half_k - P.drlevel), 28));
    std::vector<size_t> csrOff;       // prefix-sum offsets
    std::vector<uint32_t> csrPosts;   // flat posting data

    // For use64: hash map fallback (direct array too large)
    struct PostRange { size_t off; uint32_t cnt; };
    phmap::flat_hash_map<uint64_t, PostRange> postIdx64;

    csrPosts.resize(totalPostings);

    if (!use64 && HASH_SPACE > 0) {
        // Non-use64: direct array approach
        csrOff.resize(static_cast<size_t>(HASH_SPACE) + 1, 0);
        // Pass 1: count per hash
        for (int s = 0; s < NUM_SHARDS; s++)
            for (auto& [k, v] : invShards[s])
                csrOff[static_cast<uint32_t>(k) + 1] = v.size();
        // Prefix sum
        for (uint32_t h = 0; h < HASH_SPACE; h++)
            csrOff[h + 1] += csrOff[h];
        // Pass 2: scatter posting data into flat array
        for (int s = 0; s < NUM_SHARDS; s++) {
            for (auto& [k, v] : invShards[s]) {
                uint32_t h = static_cast<uint32_t>(k);
                std::memcpy(&csrPosts[csrOff[h]], v.data(),
                            v.size() * sizeof(uint32_t));
            }
            invShards[s].clear();
        }
    } else {
        // use64: hash map + flat posting array
        postIdx64.reserve(totalAfter);
        size_t pos = 0;
        for (int s = 0; s < NUM_SHARDS; s++) {
            for (auto& [k, v] : invShards[s]) {
                std::memcpy(&csrPosts[pos], v.data(),
                            v.size() * sizeof(uint32_t));
                postIdx64[k] = {pos, static_cast<uint32_t>(v.size())};
                pos += v.size();
            }
            invShards[s].clear();
        }
    }
    invShards.clear(); invShards.shrink_to_fit();
    std::cerr << "flatten CSR: " << get_sec() - tFlat << " s"
              << (use64 ? " (hash map)" : " (direct array)") << "\n";

    // Sort each genome's sketch for sequential CSR access
    double tSort = get_sec();
#pragma omp parallel for num_threads(nThreads) schedule(dynamic)
    for (int i = 0; i < N; i++) {
        if (use64) std::sort(sketches[i].h64.begin(), sketches[i].h64.end());
        else       std::sort(sketches[i].h32.begin(), sketches[i].h32.end());
    }
    std::cerr << "sort sketches: " << get_sec() - tSort << " s\n";

    // ── Phase 3: Distance via flat CSR index ─────────────────────────
    double t3 = get_sec();

    const double inv_kmer_size = 1.0 / kmerSize;
    const double p_exp   = std::exp(-static_cast<double>(kmerSize) * maxDist);
    const double minJac  = p_exp / (2.0 - p_exp);
    const double radio   = 2.0 * std::exp(maxDist * (kmerSize - 1)) - 1.0;
    std::cerr << "pruning: minJac=" << minJac << "  radio=" << radio
              << "  (mashD<" << maxDist << ", k=" << kmerSize << ")\n";

    std::vector<int> sketchSz(N);
    for (int i = 0; i < N; i++)
        sketchSz[i] = use64 ? static_cast<int>(sketches[i].h64.size())
                             : static_cast<int>(sketches[i].h32.size());

    const size_t* csrOffPtr = csrOff.data();
    const uint32_t* csrPostPtr = csrPosts.data();

    FILE* fout = fopen(outPath.c_str(), "w");
    setvbuf(fout, nullptr, _IOFBF, 1 << 24);  // 16MB kernel buffer
    fprintf(fout, "genome0\tgenome1\tcommon|size0|size1\tjaccard\tmashD\n");

#pragma omp parallel num_threads(nThreads)
    {
        std::vector<int> isect(N, 0);
        std::vector<int> stamp(N, 0);
        int ep = 0;
        std::vector<int> cand;
        cand.reserve(4096);
        std::string buf;
        buf.reserve(1 << 24);  // 16MB app buffer per thread

#pragma omp for schedule(dynamic, 64)
        for (int i = 0; i < N; i++) {
            const int s0 = sketchSz[i];
            if (__builtin_expect(s0 == 0, 0)) continue;

            cand.clear();
            ++ep;
            if (__builtin_expect(ep == INT_MAX, 0)) {
                std::memset(stamp.data(), 0, N * sizeof(int));
                ep = 1;
            }

            if (use64) {
                const auto& harr = sketches[i].h64;
                const size_t hsz = harr.size();
                for (size_t idx = 0; idx < hsz; idx++) {
                    if (__builtin_expect(idx + 1 < hsz, 1))
                        __builtin_prefetch(&harr[idx + 1], 0, 1);
                    const uint64_t hv = harr[idx];
                    auto it = postIdx64.find(hv);
                    if (__builtin_expect(it == postIdx64.end(), 0))
                        continue;
                    const uint32_t* pl = csrPostPtr + it->second.off;
                    const uint32_t plSz = it->second.cnt;
                    for (uint32_t ki = 0; ki < plSz; ki++) {
                        uint32_t j = pl[ki];
                        if (static_cast<int>(j) <= i) continue;
                        if (__builtin_expect(stamp[j] != ep, 1)) {
                            stamp[j] = ep;
                            isect[j] = 1;
                            cand.push_back(j);
                        } else {
                            isect[j]++;
                        }
                    }
                }
            } else {
                const auto& harr = sketches[i].h32;
                const size_t hsz = harr.size();
                for (size_t idx = 0; idx < hsz; idx++) {
                    if (__builtin_expect(idx + 1 < hsz, 1))
                        __builtin_prefetch(&harr[idx + 1], 0, 1);
                    const uint32_t hv = harr[idx];
                    const size_t plStart = csrOffPtr[hv];
                    const size_t plEnd   = csrOffPtr[hv + 1];
                    if (__builtin_expect(plStart == plEnd, 0)) continue;
                    for (size_t ki = plStart; ki < plEnd; ki++) {
                        uint32_t j = csrPostPtr[ki];
                        if (static_cast<int>(j) <= i) continue;
                        if (__builtin_expect(stamp[j] != ep, 1)) {
                            stamp[j] = ep;
                            isect[j] = 1;
                            cand.push_back(j);
                        } else {
                            isect[j]++;
                        }
                    }
                }
            }

            const int minCommon = std::max(1,
                static_cast<int>(std::ceil(minJac * s0)));

            for (int j : cand) {
                const int common = isect[j];
                if (common < minCommon) continue;

                const int s1 = sketchSz[j];
                const int mn = s0 < s1 ? s0 : s1;
                const int mx = s0 > s1 ? s0 : s1;
                if (__builtin_expect(mx > radio * mn, 0)) continue;

                const int denom = s0 + s1 - common;
                const double jac = static_cast<double>(common) / denom;
                const double mashD = (jac >= 1.0) ? 0.0
                    : -inv_kmer_size * std::log(2.0 * jac / (1.0 + jac));
                if (mashD < maxDist) {
                    char line[1024];
                    int n = snprintf(line, sizeof(line),
                        "%s\t%s\t%d|%d|%d\t%.6f\t%.6f\n",
                        fileList[j].c_str(), fileList[i].c_str(),
                        common, s0, s1, jac, mashD);
                    buf.append(line, n);
                }
            }

            if (buf.size() > (1 << 24)) {  // 16MB flush threshold
#pragma omp critical
                { fwrite(buf.data(), 1, buf.size(), fout); }
                buf.clear();
            }
            if (i % progress == 0)
                std::cerr << "  dist " << i << " / " << N << "\n";
        }
        if (!buf.empty()) {
#pragma omp critical
            { fwrite(buf.data(), 1, buf.size(), fout); }
        }
    }
    fclose(fout);

    double t4 = get_sec();
    std::cerr << "dist time: " << t4 - t3 << " s\n";
    std::cerr << "total: " << t4 - t0 << " s\n";

    // Skip destructors for sketches(208K vectors) + invShards(17M entries)
    // which would otherwise take 10-30s to free
    std::cerr.flush();
    std::_Exit(0);
}
