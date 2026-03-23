/**
 * ProbKMV – K Minimum Values sketch with ProbMinHash-style hashing.
 *
 * Route B: for each element x, generate a single key g(x) = Uniform(0,1)
 * from hash(x), then keep the k smallest distinct keys:
 *
 *     sketch(A) = bottom-k { g(x) : x ∈ A }
 *
 * Internal representation: a sorted double[k] array.  During warm-up
 * (size < k), new keys are inserted via binary-search + shift.  Once
 * full, most elements are rejected in O(1) by comparing against the
 * threshold (the largest value in the bottom-k).
 *
 * Jaccard estimator (standard KMV):
 *   merge two sorted sketches, take k smallest distinct values, count
 *   how many appear in both:  J ≈ common / k.
 */

#include "probmh.h"
#include "hash_int.h"

#include <immintrin.h>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>

using namespace Sketch;

// ═══════════════════════════════════════════════════════════════════════════
// DNA encoding  (same LUT as ProbMinHash4 / SetSketch)
// ═══════════════════════════════════════════════════════════════════════════

static const uint8_t ENC_LUT[256] = {
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,  0,255,  1,255,255,255,  2,255,255,255,255,255,255,255,255,
    255,255,255,255,  3,255,255,255,255,255,255,255,255,255,255,255,
    255,  0,255,  1,255,255,255,  2,255,255,255,255,255,255,255,255,
    255,255,255,255,  3,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
};
#define PMH_ENC(c)   (ENC_LUT[(uint8_t)(c)])
#define PMH_COMP(e)  ((uint8_t)((e) <= 3 ? 3 - (e) : 255))
#define PMH_VALID(e) ((e) <= 3)

static const double PMH_INV2_53 = 1.0 / (double)(1ULL << 53);

// ═══════════════════════════════════════════════════════════════════════════
// ProbKMV  –  constructor / copy / assign
// ═══════════════════════════════════════════════════════════════════════════

ProbKMV::ProbKMV(uint32_t k, int kmer_size, uint64_t seed)
    : k_(k), kmer_size_(kmer_size), seed_(seed),
      vals_(new double[k]),
      size_(0),
      threshold_(std::numeric_limits<double>::infinity())
{
    assert(k > 1);
    assert(kmer_size >= 1 && kmer_size <= 32);
    std::fill_n(vals_.get(), k, std::numeric_limits<double>::infinity());
}

ProbKMV::ProbKMV(const ProbKMV& o)
    : k_(o.k_), kmer_size_(o.kmer_size_), seed_(o.seed_),
      vals_(new double[o.k_]),
      size_(o.size_),
      threshold_(o.threshold_)
{
    std::copy(o.vals_.get(), o.vals_.get() + k_, vals_.get());
}

ProbKMV& ProbKMV::operator=(ProbKMV other) {
    std::swap(k_,         other.k_);
    std::swap(kmer_size_, other.kmer_size_);
    std::swap(seed_,      other.seed_);
    std::swap(vals_,      other.vals_);
    std::swap(size_,      other.size_);
    std::swap(threshold_, other.threshold_);
    return *this;
}

// ═══════════════════════════════════════════════════════════════════════════
// insertKey  –  maintain sorted bottom-k with deduplication
//
// Array invariant: vals_[0..size_-1] sorted ascending, distinct;
//                  vals_[size_..k_-1] = +inf.
// ═══════════════════════════════════════════════════════════════════════════

void ProbKMV::insertKey(double key) {
    if (size_ < k_) {
        double* pos = std::lower_bound(vals_.get(), vals_.get() + size_, key);
        uint32_t idx = static_cast<uint32_t>(pos - vals_.get());
        if (idx < size_ && vals_[idx] == key) return;
        std::memmove(vals_.get() + idx + 1, vals_.get() + idx,
                     (size_ - idx) * sizeof(double));
        vals_[idx] = key;
        ++size_;
        if (size_ == k_)
            threshold_ = vals_[k_ - 1];
        return;
    }

    if (key >= threshold_) return;

    double* pos = std::lower_bound(vals_.get(), vals_.get() + k_, key);
    uint32_t idx = static_cast<uint32_t>(pos - vals_.get());
    if (idx < k_ && vals_[idx] == key) return;

    std::memmove(vals_.get() + idx + 1, vals_.get() + idx,
                 (k_ - 1 - idx) * sizeof(double));
    vals_[idx] = key;
    threshold_ = vals_[k_ - 1];
}

// ═══════════════════════════════════════════════════════════════════════════
// addHash  –  hash → Uniform(0,1) key → insert into bottom-k
// ═══════════════════════════════════════════════════════════════════════════

void ProbKMV::addHash(uint64_t h) {
    uint64_t h1 = mc::murmur3_fmix(h, seed_);
    double key = (double)(h1 >> 11) * PMH_INV2_53;
    insertKey(key);
}

// ═══════════════════════════════════════════════════════════════════════════
// update  –  k-mer rolling hash with SIMD batch processing
// ═══════════════════════════════════════════════════════════════════════════

void ProbKMV::update(const char* seq, uint64_t length) {
    const int K = kmer_size_;
    if (length < static_cast<uint64_t>(K)) return;

    const uint64_t loc_seed = seed_;

    const uint64_t kmer_mask = (K == 32) ? ~0ULL : ((1ULL << (2 * K)) - 1);
    uint64_t fwd = 0, rev = 0;
    int inv = 0;
    for (int k = 0; k < K; ++k) {
        uint8_t ef = PMH_ENC(seq[k]);
        if (!PMH_VALID(ef)) inv++;
        fwd = (fwd << 2) | (PMH_VALID(ef) ? (ef & 3u) : 0u);
        uint8_t er = PMH_VALID(ef) ? PMH_COMP(ef) : 0u;
        rev = (rev >> 2) | (static_cast<uint64_t>(er & 3u) << (2 * (K - 1)));
    }

    const int      lanes   = 8;
    const uint64_t N_body  = length - static_cast<uint64_t>(K);
    const uint64_t N_batch = (N_body >= 1) ? (N_body / lanes) * lanes : 0;

    for (uint64_t i = 0; i < N_batch; i += lanes) {
        uint64_t resv[8];
        bool     lane_valid[8];

        for (int j = 0; j < lanes; ++j) {
            uint64_t pos = i + j;
            lane_valid[j] = (inv == 0);
            resv[j] = lane_valid[j] ? ((fwd <= rev) ? fwd : rev) : 0;
            uint8_t ef_out = PMH_ENC(seq[pos]);
            uint8_t ef_in  = PMH_ENC(seq[pos + K]);
            if (!PMH_VALID(ef_out)) inv--;
            if (!PMH_VALID(ef_in))  inv++;
            fwd = ((fwd << 2) | (PMH_VALID(ef_in) ? (ef_in & 3u) : 0u)) & kmer_mask;
            uint8_t er_in = PMH_VALID(ef_in) ? (PMH_COMP(ef_in) & 3u) : 0u;
            rev = (rev >> 2) | (static_cast<uint64_t>(er_in) << (2 * (K - 1)));
        }

        if (!(lane_valid[0]|lane_valid[1]|lane_valid[2]|lane_valid[3]|
              lane_valid[4]|lane_valid[5]|lane_valid[6]|lane_valid[7])) continue;

        // ── 8-lane batched murmur3_fmix (AVX-512) ───────────────────────
        uint64_t h0v[8];
#if defined(__AVX512F__) && defined(__AVX512DQ__)
        {
            __m512i vb = _mm512_loadu_si512((const void*)resv);
            __m512i vs = _mm512_set1_epi64((int64_t)loc_seed);
            __m512i va = _mm512_xor_epi64(vb, vs);
            __m512i vt = _mm512_srli_epi64(va, 33);
            vb = _mm512_xor_epi64(va, vt);
            va = _mm512_mullo_epi64(vb, _mm512_set1_epi64(0xff51afd7ed558ccdLL));
            vt = _mm512_srli_epi64(va, 33);
            vb = _mm512_xor_epi64(va, vt);
            va = _mm512_mullo_epi64(vb, _mm512_set1_epi64(0xc4ceb9fe1a85ec53LL));
            vt = _mm512_srli_epi64(va, 33);
            vb = _mm512_xor_epi64(va, vt);
            _mm512_storeu_si512(h0v, vb);
        }
#else
        for (int j = 0; j < lanes; ++j)
            h0v[j] = mc::murmur3_fmix(resv[j], loc_seed);
#endif

        // ── Per-lane: key generation + threshold filter + insert ─────────
        double cur_thresh = (size_ >= k_) ? threshold_ : 1.0;
        for (int j = 0; j < lanes; ++j) {
            if (!lane_valid[j]) continue;
#ifdef PMH_FAST_HASH
            uint64_t h1 = h0v[j];
#else
            uint64_t h1 = mc::murmur3_fmix(h0v[j], loc_seed);
#endif
            double key = (double)(h1 >> 11) * PMH_INV2_53;
            if (key >= cur_thresh) continue;
            insertKey(key);
            cur_thresh = (size_ >= k_) ? threshold_ : 1.0;
        }
    }

    // ── Remainder loop (scalar) ──────────────────────────────────────────
    for (uint64_t i = N_batch; i <= N_body; ++i) {
        if (inv == 0) {
            const uint64_t canonical = (fwd <= rev) ? fwd : rev;
#ifdef PMH_FAST_HASH
            uint64_t h1 = mc::murmur3_fmix(canonical, loc_seed);
#else
            uint64_t h0 = mc::murmur3_fmix(canonical, loc_seed);
            uint64_t h1 = mc::murmur3_fmix(h0, loc_seed);
#endif
            double key = (double)(h1 >> 11) * PMH_INV2_53;
            double thr = (size_ >= k_) ? threshold_ : 1.0;
            if (key < thr)
                insertKey(key);
        }
        if (i < N_body) {
            uint8_t ef_out = PMH_ENC(seq[i]);
            uint8_t ef_in  = PMH_ENC(seq[i + K]);
            if (!PMH_VALID(ef_out)) inv--;
            if (!PMH_VALID(ef_in))  inv++;
            fwd = ((fwd << 2) | (PMH_VALID(ef_in) ? (ef_in & 3u) : 0u)) & kmer_mask;
            uint8_t er_in = PMH_VALID(ef_in) ? (PMH_COMP(ef_in) & 3u) : 0u;
            rev = (rev >> 2) | (static_cast<uint64_t>(er_in) << (2 * (K - 1)));
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// jaccard  –  standard KMV bottom-k intersection estimator
//
// Two-pointer merge of the sorted arrays S_A and S_B.  Walk through
// the k smallest distinct values in S_A ∪ S_B, counting how many
// appear in both:  J ≈ common / k.
// ═══════════════════════════════════════════════════════════════════════════

double ProbKMV::jaccard(const ProbKMV& other) const {
    assert(k_ == other.k_);

    const double* __restrict__ a = vals_.get();
    const double* __restrict__ b = other.vals_.get();
    const uint32_t sa = size_;
    const uint32_t sb = other.size_;
    const double inf = std::numeric_limits<double>::infinity();

    uint32_t ia = 0, ib = 0;
    int distinct = 0, common = 0;

    while (static_cast<uint32_t>(distinct) < k_ && ia < sa && ib < sb) {
        if (a[ia] >= inf || b[ib] >= inf) break;
        if (a[ia] == b[ib]) {
            ++common;
            ++distinct;
            ++ia;
            ++ib;
        } else if (a[ia] < b[ib]) {
            ++distinct;
            ++ia;
        } else {
            ++distinct;
            ++ib;
        }
    }
    while (static_cast<uint32_t>(distinct) < k_ && ia < sa && a[ia] < inf) {
        ++distinct;
        ++ia;
    }
    while (static_cast<uint32_t>(distinct) < k_ && ib < sb && b[ib] < inf) {
        ++distinct;
        ++ib;
    }

    return (distinct > 0) ? static_cast<double>(common) / static_cast<double>(distinct)
                          : 0.0;
}

// ═══════════════════════════════════════════════════════════════════════════
// merge  –  sorted merge of two bottom-k lists, take k smallest distinct
// ═══════════════════════════════════════════════════════════════════════════

ProbKMV ProbKMV::merge(const ProbKMV& other) const {
    assert(k_ == other.k_ && kmer_size_ == other.kmer_size_);
    ProbKMV ret(k_, kmer_size_, seed_);

    const double* a = vals_.get();
    const double* b = other.vals_.get();
    const double inf = std::numeric_limits<double>::infinity();
    uint32_t ia = 0, ib = 0;

    while (ret.size_ < k_ && ia < size_ && ib < other.size_) {
        if (a[ia] >= inf && b[ib] >= inf) break;
        double next;
        if (a[ia] == b[ib]) {
            next = a[ia];
            ++ia;
            ++ib;
        } else if (a[ia] < b[ib]) {
            next = a[ia++];
        } else {
            next = b[ib++];
        }
        ret.vals_[ret.size_++] = next;
    }
    while (ret.size_ < k_ && ia < size_ && a[ia] < inf) {
        ret.vals_[ret.size_++] = a[ia++];
    }
    while (ret.size_ < k_ && ib < other.size_ && b[ib] < inf) {
        ret.vals_[ret.size_++] = b[ib++];
    }

    ret.threshold_ = (ret.size_ >= k_) ? ret.vals_[k_ - 1] : inf;
    return ret;
}

// ═══════════════════════════════════════════════════════════════════════════
// printSketch
// ═══════════════════════════════════════════════════════════════════════════

void ProbKMV::printSketch() const {
    std::fprintf(stdout,
                 "ProbKMV k=%u kmer=%d fill=%u/%u vals[0..19]: ",
                 k_, kmer_size_, size_, k_);
    for (uint32_t i = 0; i < k_ && i < 20; ++i)
        std::fprintf(stdout, "%.4g ", vals_[i]);
    if (k_ > 20) std::fprintf(stdout, "...");
    std::fprintf(stdout, "\n");
}

#undef PMH_ENC
#undef PMH_COMP
#undef PMH_VALID
