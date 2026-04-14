/**
 * InvertedIndex.h – generic inverted-index module for sketch-based all-to-all
 * distance computation.
 *
 * Header-only, templated on key type (uint32_t or uint64_t).
 *
 * Workflow:
 *   1. Caller builds per-thread local inverted indices during sketch construction.
 *   2. buildCSRIndex()  – 64-shard parallel merge + singleton removal + CSR flatten.
 *   3. computeDistances() – posting-list traversal with stamp/epoch counting,
 *      caller-provided Jaccard lambda, minCommon pruning, and Mash distance output.
 */

#ifndef _INVERTED_INDEX_H_
#define _INVERTED_INDEX_H_

#include "phmap.h"
#include "common.h"

#include <omp.h>
#include <cstring>
#include <cmath>
#include <climits>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <sys/stat.h>

namespace Sketch {

// ─────────────────────────────────────────────────────────────────────────────
// Data structure holding the CSR-packed inverted index.
// ─────────────────────────────────────────────────────────────────────────────
template<typename KeyT>
struct InvertedIndex {
    struct PostRange { size_t off; uint32_t cnt; };

    phmap::flat_hash_map<KeyT, PostRange> postIdx;
    std::vector<uint32_t>                 csrPosts;
    size_t                                totalPostings = 0;
};

// ─────────────────────────────────────────────────────────────────────────────
// Phase 2: merge thread-local maps → singleton removal → CSR flatten
// ─────────────────────────────────────────────────────────────────────────────
template<typename KeyT>
InvertedIndex<KeyT> buildCSRIndex(
    std::vector<phmap::flat_hash_map<KeyT, std::vector<uint32_t>>>& threadIdx,
    int nThreads)
{
    constexpr int NUM_SHARDS = 64;
    const KeyT SHARD_MASK = static_cast<KeyT>(NUM_SHARDS - 1);
    const int tCount = static_cast<int>(threadIdx.size());

    std::vector<phmap::flat_hash_map<KeyT, std::vector<uint32_t>>> invShards(NUM_SHARDS);

    // Phase 2.1: thread-local re-shard (no locks)
    {
        std::vector<std::vector<phmap::flat_hash_map<KeyT, std::vector<uint32_t>>>> localShards(
            tCount, std::vector<phmap::flat_hash_map<KeyT, std::vector<uint32_t>>>(NUM_SHARDS));

        #pragma omp parallel for num_threads(nThreads) schedule(static)
        for (int tid = 0; tid < tCount; ++tid)
        {
            for (auto& [key, vec] : threadIdx[tid]) {
                int shard = static_cast<int>(key & SHARD_MASK);
                auto& shardMap = localShards[tid][shard];
                auto [it, ins] = shardMap.try_emplace(key, std::move(vec));
                if (!ins) {
                    auto& dst = it->second;
                    dst.reserve(dst.size() + vec.size());
                    dst.insert(dst.end(), vec.begin(), vec.end());
                }
            }
            phmap::flat_hash_map<KeyT, std::vector<uint32_t>>().swap(threadIdx[tid]);
        }
        threadIdx.clear();

        // Phase 2.2: shard merge across threads
        #pragma omp parallel for num_threads(nThreads) schedule(dynamic, 1)
        for (int s = 0; s < NUM_SHARDS; s++) {
            auto& merged = invShards[s];
            for (int tid = 0; tid < tCount; ++tid) {
                auto& src = localShards[tid][s];
                for (auto& [key, vec] : src) {
                    auto [it, ins] = merged.try_emplace(key, std::move(vec));
                    if (!ins) {
                        auto& dst = it->second;
                        dst.reserve(dst.size() + vec.size());
                        dst.insert(dst.end(), vec.begin(), vec.end());
                    }
                }
                phmap::flat_hash_map<KeyT, std::vector<uint32_t>>().swap(src);
            }
        }
    }

    size_t totalUnique = 0;
    for (auto& sh : invShards) totalUnique += sh.size();
    std::cerr << "merge index: " << totalUnique << " unique keys" << std::endl;

    // Singleton removal + posting count in one pass
    size_t totalAfter = 0;
    size_t totalPostings = 0;
    #pragma omp parallel for num_threads(nThreads) reduction(+:totalAfter,totalPostings)
    for (int s = 0; s < NUM_SHARDS; s++) {
        for (auto it = invShards[s].begin(); it != invShards[s].end(); ) {
            if (it->second.size() <= 1) it = invShards[s].erase(it);
            else {
                totalPostings += it->second.size();
                ++it;
                totalAfter++;
            }
        }
    }
    std::cerr << "singleton removal: " << totalUnique << " -> " << totalAfter << std::endl;

    // Flatten into CSR
    InvertedIndex<KeyT> idx;
    idx.totalPostings = totalPostings;
    idx.postIdx.reserve(totalAfter);
    idx.csrPosts.resize(totalPostings);

    size_t cursor = 0;
    for (int s = 0; s < NUM_SHARDS; s++) {
        for (auto& [key, vec] : invShards[s]) {
            idx.postIdx[key] = {cursor, static_cast<uint32_t>(vec.size())};
            std::memcpy(&idx.csrPosts[cursor], vec.data(),
                        vec.size() * sizeof(uint32_t));
            cursor += vec.size();
        }
        phmap::flat_hash_map<KeyT, std::vector<uint32_t>>().swap(invShards[s]);
    }

    std::cerr << "CSR flatten: " << totalPostings << " postings" << std::endl;
    return idx;
}

// ─────────────────────────────────────────────────────────────────────────────
// Phase 3: all-to-all distance via CSR posting traversal.
//
// JaccardFn:    double fn(int common, int s_i, int s_j)  → Jaccard similarity
// MinCommonFn:  int fn(int s_i) → minimum common count to pass distance filter
//
// FastKMV/Kssd call with:
//   jaccardFn   = [](int c, int s0, int s1){ return (double)c/(s0+s1-c); }
//   minCommonFn = [minJac](int s0){ return max(1,(int)ceil(minJac*s0)); }
//
// BinDash call with:
//   jaccardFn   = [p,n](int c,int,int){ return (c/n-p)/(1-p); }
//   minCommonFn = [minSame](int){ return minSame; }   // global constant
// ─────────────────────────────────────────────────────────────────────────────
template<typename KeyT, typename JaccardFn, typename MinCommonFn>
void computeDistances(
    const InvertedIndex<KeyT>&              idx,
    const std::vector<std::vector<KeyT>>&   skKeys,
    const std::vector<int>&                 sketchSizes,
    const std::vector<std::string>&         fileList,
    int                                     N,
    int                                     kmerSize,
    double                                  maxDist,
    JaccardFn                               jaccardFn,
    MinCommonFn                             minCommonFn,
    const std::string&                      outputPath,
    int                                     nThreads)
{
    const bool useMashDist = (kmerSize > 0);
    const uint32_t* csrPtr = idx.csrPosts.data();
    double minJacByDist = 1.0 - maxDist;
    if (useMashDist) {
        const double p_exp = std::exp(-static_cast<double>(kmerSize) * maxDist);
        minJacByDist = p_exp / (2.0 - p_exp);
    }

    // Size-ratio pruning (only meaningful for Mash distance mode)
    const double radio = useMashDist
        ? 2.0 * std::exp(maxDist * (kmerSize - 1)) - 1.0
        : 1e18;

    // If outputPath is a directory, append a default filename
    std::string finalPath = outputPath;
    struct stat st;
    if (stat(finalPath.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
        if (finalPath.back() != '/') finalPath += '/';
        finalPath += "res.dist";
    }

    FILE* fout = fopen(finalPath.c_str(), "w");
    if (!fout) {
        std::cerr << "ERROR: cannot open output file: " << finalPath << std::endl;
        return;
    }
    setvbuf(fout, nullptr, _IOFBF, 1 << 24);
    std::cerr << "output: " << finalPath << std::endl;

    int progress = N / 20;
    if (progress < 1) progress = 1;

    #pragma omp parallel num_threads(nThreads)
    {
        std::vector<int> isect(N, 0);
        std::vector<int> stamp(N, 0);
        int ep = 0;
        std::vector<int> cand;
        cand.reserve(4096);
        std::string buf;
        buf.reserve(1 << 24);

        #pragma omp for schedule(dynamic, 64)
        for (int i = 0; i < N; i++) {
            const int s0 = sketchSizes[i];
            if (__builtin_expect(s0 == 0, 0)) continue;

            const int minCommon = minCommonFn(s0);

            cand.clear();
            ++ep;
            if (__builtin_expect(ep == INT_MAX, 0)) {
                std::memset(stamp.data(), 0, N * sizeof(int));
                ep = 1;
            }

            const auto& keys = skKeys[i];
            const size_t ksz = keys.size();
            for (size_t ki = 0; ki < ksz; ki++) {
                if (__builtin_expect(ki + 1 < ksz, 1))
                    __builtin_prefetch(&keys[ki + 1], 0, 1);
                auto it = idx.postIdx.find(keys[ki]);
                if (__builtin_expect(it == idx.postIdx.end(), 0)) continue;
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

            for (int j : cand) {
                const int common = isect[j];
                if (common < minCommon) continue;

                const int s1 = sketchSizes[j];
                const int mn = s0 < s1 ? s0 : s1;
                const int mx = s0 > s1 ? s0 : s1;
                if (__builtin_expect(mx > radio * mn, 0)) continue;

                double jac = jaccardFn(common, s0, s1);
                if (jac < 0.0) jac = 0.0;
                if (jac > 1.0) jac = 1.0;

                if (jac > minJacByDist) {
                    double dist = 1.0 - jac;
                    if (useMashDist) {
                        dist = (jac >= 1.0) ? 0.0
                            : -std::log(2.0 * jac / (1.0 + jac))
                              / static_cast<double>(kmerSize);
                    }
                    char line[1024];
                    int n = snprintf(line, sizeof(line), "%s\t%s\t%.6f\n",
                        fileList[i].c_str(), fileList[j].c_str(), dist);
                    buf.append(line, static_cast<size_t>(n));
                }
            }

            if (buf.size() > (1 << 24)) {
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
}

} // namespace Sketch

#endif // _INVERTED_INDEX_H_
