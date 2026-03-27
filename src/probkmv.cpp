/**
 * ProbKMV – K Minimum Values sketch with ntHash rolling + fmix finalizer.
 *
 * For each k-mer x, produce a canonical hash via ntHash
 * (Mohamadi et al., Bioinformatics 2016), then apply one round of
 * murmur3 fmix as a finalizer.  The 53-bit key  g(x) = fmix(H) >> 11
 * is kept in a sorted uint64_t[k] bottom-k structure.
 *
 * Rolling hash update per position (3-cycle critical path):
 *   H_fwd' = rol(H_fwd, 1) ⊕ rol(seed[out], k) ⊕ seed[in]
 *   H_rc'  = rol(H_rc,  1) ⊕ rol(seedC[out], k) ⊕ seedC[in]
 *   canonical = H_fwd ⊕ H_rc   (branch-free, no min/cmov)
 *
 * Entire pipeline stays in the integer domain; SIMD 8-wide fmix
 * finalizer + VPCMPUQ threshold filter.
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
#include <climits>

using namespace Sketch;

// ═══════════════════════════════════════════════════════════════════════════
// DNA encoding via bit operations (no 256-byte LUT)
//
//   (c >> 1) & 3:  A→0, C→1, T→2, G→3  (case-insensitive)
//   complement:    e ^ 2   (A↔T = 0↔2, C↔G = 1↔3)
//
// Validity is checked separately; the pre-encoding pass marks non-ACGT
// bytes as 0xFF so the rolling loop can test (e <= 3) as before.
// ═══════════════════════════════════════════════════════════════════════════

static inline uint8_t pmh_enc_or_invalid(uint8_t c) {
    const uint8_t cu = c & 0xDF;
    if (cu == 'A' || cu == 'C' || cu == 'G' || cu == 'T')
        return (c >> 1) & 3;
    return 0xFF;
}
#define PMH_VALID(e) ((e) <= 3)

// ═══════════════════════════════════════════════════════════════════════════
// ntHash seed tables  (values from Mohamadi et al., Bioinformatics 2016)
//
// Correct ntHash rolling formulas (derived from paper):
//
//   H_fwd = XOR_{k=0}^{K-1} rol(seed_fwd[enc[k]], K-1-k)
//   H_rc  = XOR_{k=0}^{K-1} rol(seed_rc[enc[k]],  k)
//
//   Forward rolling:  H_fwd' = rol(H_fwd, 1) ⊕ rol(seed_fwd[out], K) ⊕ seed_fwd[in]
//   RC rolling:       H_rc'  = ror(H_rc,  1) ⊕ ror(seed_rc[out],  1) ⊕ rol(seed_rc[in], K-1)
//
//   Canonical = min(H_fwd, H_rc)   [NOT xor — xor maps all palindromic k-mers to 0]
// ═══════════════════════════════════════════════════════════════════════════

static const uint64_t NT_SEED_FWD[4] = {
    0x3c8bfbb395c60474ULL,   // A = 0
    0x3193c18562a02b4cULL,   // C = 1
    0x295549f54be24456ULL,   // T = 2
    0x20323ed082572324ULL,   // G = 3
};
static const uint64_t NT_SEED_RC[4] = {
    NT_SEED_FWD[2],          // comp(A=0) → T seed
    NT_SEED_FWD[3],          // comp(C=1) → G seed
    NT_SEED_FWD[0],          // comp(T=2) → A seed
    NT_SEED_FWD[1],          // comp(G=3) → C seed
};

static inline uint64_t rol64(uint64_t x, unsigned n) {
    n &= 63;
    return (x << n) | (x >> ((64 - n) & 63));
}

// ═══════════════════════════════════════════════════════════════════════════
// ProbKMV  –  constructor / copy / assign
// ═══════════════════════════════════════════════════════════════════════════

ProbKMV::ProbKMV(uint32_t k, int kmer_size, uint64_t seed)
    : k_(k), kmer_size_(kmer_size), seed_(seed),
      vals_(new uint64_t[k]),
      size_(0),
      threshold_(UINT64_MAX)
{
    assert(k > 1);
    assert(kmer_size >= 1 && kmer_size <= 32);
    std::fill_n(vals_.get(), k, UINT64_MAX);
}

ProbKMV::ProbKMV(const ProbKMV& o)
    : k_(o.k_), kmer_size_(o.kmer_size_), seed_(o.seed_),
      vals_(new uint64_t[o.k_]),
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
// insertKey  –  maintain sorted bottom-k (uint64_t) with deduplication
//
// Array invariant: vals_[0..size_-1] sorted ascending, distinct;
//                  vals_[size_..k_-1] = UINT64_MAX.
// ═══════════════════════════════════════════════════════════════════════════

void ProbKMV::insertKey(uint64_t key) {
    if (size_ < k_) {
        uint64_t* pos = std::lower_bound(vals_.get(), vals_.get() + size_, key);
        uint32_t idx = static_cast<uint32_t>(pos - vals_.get());
        if (idx < size_ && vals_[idx] == key) return;
        std::memmove(vals_.get() + idx + 1, vals_.get() + idx,
                     (size_ - idx) * sizeof(uint64_t));
        vals_[idx] = key;
        ++size_;
        if (size_ == k_)
            threshold_ = vals_[k_ - 1];
        return;
    }

    if (key >= threshold_) return;

    uint64_t* pos = std::lower_bound(vals_.get(), vals_.get() + k_, key);
    uint32_t idx = static_cast<uint32_t>(pos - vals_.get());
    if (idx < k_ && vals_[idx] == key) return;

    std::memmove(vals_.get() + idx + 1, vals_.get() + idx,
                 (k_ - 1 - idx) * sizeof(uint64_t));
    vals_[idx] = key;
    threshold_ = vals_[k_ - 1];
}

// ═══════════════════════════════════════════════════════════════════════════
// addHash  –  hash → 53-bit integer key → insert into bottom-k
// ═══════════════════════════════════════════════════════════════════════════

void ProbKMV::addHash(uint64_t h) {
    uint64_t h1 = mc::murmur3_fmix(h, seed_);
    insertKey(h1 >> 11);
}

// ═══════════════════════════════════════════════════════════════════════════
// update  –  ntHash rolling + SIMD fmix finalizer + threshold filter
//
// Phase 0 – Bulk nucleotide encoding (VPSHUFB, 64 bases/cycle).
//
// Phase 1 – ntHash canonical rolling hash:
//   H_fwd and H_rc are maintained by two INDEPENDENT rolling chains,
//   allowing out-of-order CPUs to execute them in parallel (~4 cycles
//   per k-mer vs ~10 for the old 2-bit-shift approach).
//
//   canonical = min(H_fwd, H_rc)  — standard ntHash canonical, avoids
//   the zero-collision problem of XOR on palindromic k-mers.
//
//   Two rounds of murmur3 fmix as finalizer (SIMD 8-wide) for full
//   avalanche; -DPMH_FAST_HASH reduces to one round.
// ═══════════════════════════════════════════════════════════════════════════

void ProbKMV::update(const char* seq, uint64_t length) {
    const int K = kmer_size_;
    if (length < static_cast<uint64_t>(K)) return;

    const uint64_t loc_seed = seed_;

    // ── Phase 0: SIMD bulk sequence encoding ────────────────────────────
    uint8_t* enc;
    std::unique_ptr<uint8_t[]> enc_storage(new uint8_t[length]);
    enc = enc_storage.get();

    uint64_t p = 0;
#if defined(__AVX512BW__)
    {
        const __m512i vlut = _mm512_broadcast_i32x4(_mm_setr_epi8(
            -1,  0, -1,  1,  2, -1, -1,  3, -1, -1, -1, -1, -1, -1, -1, -1));
        const __m512i vmask_lo = _mm512_set1_epi8(0x0F);
        for (; p + 64 <= length; p += 64) {
            __m512i vraw = _mm512_loadu_si512(seq + p);
            __m512i vnib = _mm512_and_si512(vraw, vmask_lo);
            __m512i venc = _mm512_shuffle_epi8(vlut, vnib);
            _mm512_storeu_si512(enc + p, venc);
        }
    }
#elif defined(__AVX2__)
    {
        const __m256i vlut = _mm256_broadcastsi128_si256(_mm_setr_epi8(
            -1,  0, -1,  1,  2, -1, -1,  3, -1, -1, -1, -1, -1, -1, -1, -1));
        const __m256i vmask_lo = _mm256_set1_epi8(0x0F);
        for (; p + 32 <= length; p += 32) {
            __m256i vraw = _mm256_loadu_si256((const __m256i*)(seq + p));
            __m256i vnib = _mm256_and_si256(vraw, vmask_lo);
            __m256i venc = _mm256_shuffle_epi8(vlut, vnib);
            _mm256_storeu_si256((__m256i*)(enc + p), venc);
        }
    }
#endif
    for (; p < length; ++p)
        enc[p] = pmh_enc_or_invalid((uint8_t)seq[p]);

    // ── Phase 1: ntHash rolling hash ────────────────────────────────────

    // Precompute rolling constants:
    //   sf_k[i]    = rol(seed_fwd[i], K)   — fwd remove (rol by K)
    //   sc_ror1[i] = ror(seed_rc[i],  1)   — rc  remove (ror by 1 = rol by 63)
    //   sc_km1[i]  = rol(seed_rc[i],  K-1) — rc  add    (rol by K-1)
    uint64_t sf_k[4], sc_ror1[4], sc_km1[4];
    const unsigned km1 = static_cast<unsigned>(K > 1 ? K - 1 : 63);
    for (int i = 0; i < 4; ++i) {
        sf_k[i]    = rol64(NT_SEED_FWD[i], static_cast<unsigned>(K));
        sc_ror1[i] = rol64(NT_SEED_RC[i],  63u);   // ror(x,1) = rol(x,63)
        sc_km1[i]  = rol64(NT_SEED_RC[i],  km1);
    }

    // Init: H_fwd = XOR_{k=0}^{K-1} rol(seed_fwd[enc[k]], K-1-k)
    //        H_rc  = XOR_{k=0}^{K-1} rol(seed_rc[enc[k]],  k)
    uint64_t fwd_h = 0, rc_h = 0;
    int inv = 0;
    for (int k = 0; k < K; ++k) {
        uint8_t e = enc[k];
        if (e > 3) { ++inv; e = 0; }
        fwd_h ^= rol64(NT_SEED_FWD[e], static_cast<unsigned>(K - 1 - k));
        rc_h  ^= rol64(NT_SEED_RC[e],  static_cast<unsigned>(k));      // FIX: k, not K-1-k
    }

    const int      lanes   = 8;
    const uint64_t N_body  = length - static_cast<uint64_t>(K);
    const uint64_t N_batch = (N_body >= 1) ? (N_body / lanes) * lanes : 0;

#if defined(__AVX512F__) && defined(__AVX512DQ__)
    const __m512i vs  = _mm512_set1_epi64((int64_t)loc_seed);
    const __m512i vc1 = _mm512_set1_epi64(0xff51afd7ed558ccdLL);
    const __m512i vc2 = _mm512_set1_epi64(0xc4ceb9fe1a85ec53LL);
    __m512i vthresh = _mm512_set1_epi64((int64_t)threshold_);
#endif

    for (uint64_t i = 0; i < N_batch; i += lanes) {
        uint64_t resv[8];
        bool     lane_valid[8];

        for (int j = 0; j < lanes; ++j) {
            const uint64_t pos = i + j;

            lane_valid[j] = (inv == 0);
            // canonical = min(fwd_h, rc_h) — NOT xor (xor→0 for palindromes)
            resv[j] = lane_valid[j] ? (fwd_h < rc_h ? fwd_h : rc_h) : 0;

            const uint8_t e_out = enc[pos];
            const uint8_t e_in  = enc[pos + K];
            if (e_out > 3) --inv;
            if (e_in  > 3) ++inv;
            const uint8_t oi = (e_out <= 3) ? e_out : 0;
            const uint8_t ii = (e_in  <= 3) ? e_in  : 0;

            // fwd: rol(H,1) ^ rol(seed_fwd[out], K) ^ seed_fwd[in]
            fwd_h = rol64(fwd_h, 1)  ^ sf_k[oi]   ^ NT_SEED_FWD[ii];
            // rc:  ror(H,1) ^ ror(seed_rc[out], 1) ^ rol(seed_rc[in], K-1)
            rc_h  = rol64(rc_h, 63u) ^ sc_ror1[oi] ^ sc_km1[ii];     // FIX
        }

#if defined(__AVX512F__) && defined(__AVX512DQ__)
        __mmask8 valid_mask = 0;
        for (int j = 0; j < lanes; ++j)
            if (lane_valid[j]) valid_mask |= (1u << j);
        if (!valid_mask) continue;

        // fmix finalizer round 1
        __m512i vb = _mm512_loadu_si512((const void*)resv);
        __m512i va = _mm512_xor_epi64(vb, vs);
        __m512i vt = _mm512_srli_epi64(va, 33);
        vb = _mm512_xor_epi64(va, vt);
        va = _mm512_mullo_epi64(vb, vc1);
        vt = _mm512_srli_epi64(va, 33);
        vb = _mm512_xor_epi64(va, vt);
        va = _mm512_mullo_epi64(vb, vc2);
        vt = _mm512_srli_epi64(va, 33);
        vb = _mm512_xor_epi64(va, vt);

#ifndef PMH_FAST_HASH
        // fmix finalizer round 2 (full avalanche for LSH banding)
        va = _mm512_xor_epi64(vb, vs);
        vt = _mm512_srli_epi64(va, 33);
        vb = _mm512_xor_epi64(va, vt);
        va = _mm512_mullo_epi64(vb, vc1);
        vt = _mm512_srli_epi64(va, 33);
        vb = _mm512_xor_epi64(va, vt);
        va = _mm512_mullo_epi64(vb, vc2);
        vt = _mm512_srli_epi64(va, 33);
        vb = _mm512_xor_epi64(va, vt);
#endif

        __m512i vkeys = _mm512_srli_epi64(vb, 11);
        __mmask8 pass_mask = _mm512_mask_cmp_epu64_mask(
            valid_mask, vkeys, vthresh, _MM_CMPINT_LT);

        if (pass_mask) {
            uint64_t keyv[8];
            _mm512_storeu_si512(keyv, vkeys);
            while (pass_mask) {
                int j = __builtin_ctz(pass_mask);
                pass_mask &= pass_mask - 1;
                insertKey(keyv[j]);
            }
            vthresh = _mm512_set1_epi64((int64_t)threshold_);
        }
#else
        // ── Scalar fallback ─────────────────────────────────────────────
        if (!(lane_valid[0]|lane_valid[1]|lane_valid[2]|lane_valid[3]|
              lane_valid[4]|lane_valid[5]|lane_valid[6]|lane_valid[7])) continue;

        for (int j = 0; j < lanes; ++j) {
            if (!lane_valid[j]) continue;
            uint64_t h0 = mc::murmur3_fmix(resv[j], loc_seed);
#ifdef PMH_FAST_HASH
            uint64_t h1 = h0;
#else
            uint64_t h1 = mc::murmur3_fmix(h0, loc_seed);
#endif
            uint64_t key = h1 >> 11;
            if (key >= threshold_) continue;
            insertKey(key);
        }
#endif
    }

    // ── Remainder loop ──────────────────────────────────────────────────
    for (uint64_t i = N_batch; i <= N_body; ++i) {
        if (inv == 0) {
            uint64_t canon = (fwd_h < rc_h) ? fwd_h : rc_h;   // FIX: min not xor
            uint64_t h0 = mc::murmur3_fmix(canon, loc_seed);
#ifdef PMH_FAST_HASH
            uint64_t h1 = h0;
#else
            uint64_t h1 = mc::murmur3_fmix(h0, loc_seed);
#endif
            uint64_t key = h1 >> 11;
            if (key < threshold_)
                insertKey(key);
        }
        if (i < N_body) {
            const uint8_t e_out = enc[i];
            const uint8_t e_in  = enc[i + K];
            if (e_out > 3) --inv;
            if (e_in  > 3) ++inv;
            const uint8_t oi = (e_out <= 3) ? e_out : 0;
            const uint8_t ii = (e_in  <= 3) ? e_in  : 0;
            fwd_h = rol64(fwd_h, 1)  ^ sf_k[oi]   ^ NT_SEED_FWD[ii];
            rc_h  = rol64(rc_h, 63u) ^ sc_ror1[oi] ^ sc_km1[ii];    // FIX
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

    const uint64_t* __restrict__ a = vals_.get();
    const uint64_t* __restrict__ b = other.vals_.get();
    const uint32_t sa = size_;
    const uint32_t sb = other.size_;

    uint32_t ia = 0, ib = 0;
    int distinct = 0, common = 0;

    while (static_cast<uint32_t>(distinct) < k_ && ia < sa && ib < sb) {
        if (a[ia] == UINT64_MAX || b[ib] == UINT64_MAX) break;
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
    while (static_cast<uint32_t>(distinct) < k_ && ia < sa && a[ia] != UINT64_MAX) {
        ++distinct;
        ++ia;
    }
    while (static_cast<uint32_t>(distinct) < k_ && ib < sb && b[ib] != UINT64_MAX) {
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

    const uint64_t* a = vals_.get();
    const uint64_t* b = other.vals_.get();
    uint32_t ia = 0, ib = 0;

    while (ret.size_ < k_ && ia < size_ && ib < other.size_) {
        if (a[ia] == UINT64_MAX && b[ib] == UINT64_MAX) break;
        uint64_t next;
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
    while (ret.size_ < k_ && ia < size_ && a[ia] != UINT64_MAX) {
        ret.vals_[ret.size_++] = a[ia++];
    }
    while (ret.size_ < k_ && ib < other.size_ && b[ib] != UINT64_MAX) {
        ret.vals_[ret.size_++] = b[ib++];
    }

    ret.threshold_ = (ret.size_ >= k_) ? ret.vals_[k_ - 1] : UINT64_MAX;
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
        std::fprintf(stdout, "%016lx ", (unsigned long)vals_[i]);
    if (k_ > 20) std::fprintf(stdout, "...");
    std::fprintf(stdout, "\n");
}

#undef PMH_VALID
