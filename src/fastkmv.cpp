/**
 * FastKMV – K Minimum Values sketch with ntHash rolling + fmix finalizer.
 *
 * Lazy-sort warmup: the first 2k keys are appended unsorted in O(1).
 * When the buffer fills, compactify() sorts, deduplicates, and truncates
 * to the bottom-k in one pass, producing a sorted array with exact
 * threshold — identical to the original design from that point on.
 *
 * Warmup (sorted_==false): O(1) append into 2k buffer; at 2k entries,
 * compactify() sorts, dedupes, truncates to bottom-k (may repeat if still <k).
 * Steady state: sorted bottom-k, O(log k) lower_bound + O(k) memmove;
 * threshold_ == v[k-1] is exact.
 *
 * Rolling hash: ntHash (Mohamadi et al. 2016), canonical = min(fwd, rc).
 * Finalizer: murmur3 fmix (default 1 round; compile with -DFASTKMV_DOUBLE_FMUX
 * for 2 rounds), 8-wide AVX-512 when available.
 *
 * Jaccard estimator (standard KMV):
 *   merge two sorted sketches, take k smallest distinct values, count
 *   how many appear in both:  J ≈ common / distinct.
 */

#include "fastkmv.h"
#include "hash_int.h"

// Default: single fmix (faster). Define FASTKMV_DOUBLE_FMUX for two rounds (legacy).
#ifdef FASTKMV_DOUBLE_FMUX
#define FASTKMV_FMUX_ROUNDS 2
#else
#define FASTKMV_FMUX_ROUNDS 1
#endif

#include <immintrin.h>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <cassert>
#include <climits>
#include <memory>

using namespace Sketch;

// ═══════════════════════════════════════════════════════════════════════════
// DNA encoding via bit operations (no 256-byte LUT)
// ═══════════════════════════════════════════════════════════════════════════

static inline uint8_t fkmv_enc_or_invalid(uint8_t c) {
    const uint8_t cu = c & 0xDF;
    if (cu == 'A' || cu == 'C' || cu == 'G' || cu == 'T')
        return (c >> 1) & 3;
    return 0xFF;
}
#define FKMV_VALID(e) ((e) <= 3)

// ═══════════════════════════════════════════════════════════════════════════
// ntHash seed tables  (Mohamadi et al., Bioinformatics 2016)
// ═══════════════════════════════════════════════════════════════════════════

static const uint64_t NT_SEED_FWD[4] = {
    0x3c8bfbb395c60474ULL,
    0x3193c18562a02b4cULL,
    0x295549f54be24456ULL,
    0x20323ed082572324ULL,
};
static const uint64_t NT_SEED_RC[4] = {
    NT_SEED_FWD[2], NT_SEED_FWD[3], NT_SEED_FWD[0], NT_SEED_FWD[1],
};

static inline uint64_t rol64(uint64_t x, unsigned n) {
    n &= 63;
    return (x << n) | (x >> ((64 - n) & 63));
}

#ifdef __AVX2__
static inline __m256i avx2_mullo_epi64(__m256i a, __m256i b) {
    __m256i a_hi  = _mm256_srli_epi64(a, 32);
    __m256i b_hi  = _mm256_srli_epi64(b, 32);
    __m256i lo_lo = _mm256_mul_epu32(a, b);
    __m256i cross = _mm256_add_epi64(_mm256_mul_epu32(a_hi, b),
                                      _mm256_mul_epu32(a, b_hi));
    return _mm256_add_epi64(lo_lo, _mm256_slli_epi64(cross, 32));
}
#endif

// ═══════════════════════════════════════════════════════════════════════════
// FastKMV  –  constructor / copy / assign
// ═══════════════════════════════════════════════════════════════════════════

FastKMV::FastKMV(uint32_t k, int kmer_size, uint64_t seed)
    : k_(k), kmer_size_(kmer_size), seed_(seed),
      buf_cap_(k * 2),
      vals_(new uint64_t[k * 2]),
      size_(0),
      threshold_(UINT64_MAX),
      sorted_(false)
{
    assert(k > 1);
    assert(kmer_size >= 1 && kmer_size <= 32);
    // Warmup: append-only until buf full → compactify (see insertKey).
    // Sentinels UINT64_MAX are applied in compactify / merge, not here.
}

FastKMV::FastKMV(const FastKMV& o)
    : k_(o.k_), kmer_size_(o.kmer_size_), seed_(o.seed_),
      buf_cap_(o.buf_cap_),
      vals_(new uint64_t[o.buf_cap_]),
      size_(o.size_),
      threshold_(o.threshold_),
      sorted_(o.sorted_)
{
    std::copy(o.vals_.get(), o.vals_.get() + size_, vals_.get());
    if (sorted_ && size_ < k_)
        std::fill(vals_.get() + size_, vals_.get() + k_, UINT64_MAX);
}

FastKMV& FastKMV::operator=(FastKMV other) {
    std::swap(k_,         other.k_);
    std::swap(kmer_size_, other.kmer_size_);
    std::swap(seed_,      other.seed_);
    std::swap(buf_cap_,   other.buf_cap_);
    std::swap(vals_,      other.vals_);
    std::swap(size_,      other.size_);
    std::swap(threshold_, other.threshold_);
    std::swap(sorted_,    other.sorted_);
    return *this;
}

// ═══════════════════════════════════════════════════════════════════════════
// insertKey
//
//   Phase A — warmup (!sorted_): O(1) append; threshold_ is UINT64_MAX.
//   Phase B — sorted but unfilled (sorted_ && size_<k): lower_bound + memmove.
//   Phase C — full KMV (sorted_ && size_==k): threshold + lower_bound + memmove.
// ═══════════════════════════════════════════════════════════════════════════

void FastKMV::insertKey(uint64_t key) {
    if (key >= threshold_) return;

    if (!sorted_) {
        vals_[size_++] = key;
        if (size_ >= buf_cap_)
            compactify();
        return;
    }

    uint64_t* v = vals_.get();

    if (size_ < k_) {
        uint64_t* pos = std::lower_bound(v, v + size_, key);
        uint32_t idx = static_cast<uint32_t>(pos - v);
        if (idx < size_ && v[idx] == key) return;
        std::memmove(v + idx + 1, v + idx, (size_ - idx) * sizeof(uint64_t));
        v[idx] = key;
        ++size_;
        if (size_ == k_)
            threshold_ = v[k_ - 1];
        return;
    }

    uint64_t* pos = std::lower_bound(v, v + k_, key);
    uint32_t idx = static_cast<uint32_t>(pos - v);
    if (idx < k_ && v[idx] == key) return;
    std::memmove(v + idx + 1, v + idx, (k_ - 1 - idx) * sizeof(uint64_t));
    v[idx] = key;
    threshold_ = v[k_ - 1];
}

// ═══════════════════════════════════════════════════════════════════════════
// compactify  –  sort, deduplicate, truncate to k
// ═══════════════════════════════════════════════════════════════════════════

void FastKMV::compactify() const {
    uint64_t* v = vals_.get();
    std::sort(v, v + size_);

    uint32_t w = 0;
    for (uint32_t r = 0; r < size_; ++r) {
        if (r == 0 || v[r] != v[r - 1]) {
            v[w++] = v[r];
            if (w == k_) break;
        }
    }
    size_ = w;
    threshold_ = (size_ >= k_) ? v[k_ - 1] : UINT64_MAX;
    std::fill(v + size_, v + k_, UINT64_MAX);
    sorted_ = true;
}

void FastKMV::ensureSorted() const {
    if (!sorted_)
        compactify();
}

// ═══════════════════════════════════════════════════════════════════════════
// addHash  –  hash → 53-bit integer key → insert into bottom-k
// ═══════════════════════════════════════════════════════════════════════════

void FastKMV::addHash(uint64_t h) {
    uint64_t h1 = mc::murmur3_fmix(h, seed_);
    insertKey(h1 >> 11);
}

// ═══════════════════════════════════════════════════════════════════════════
// update  –  ntHash rolling + SIMD fmix finalizer + threshold filter
// ═══════════════════════════════════════════════════════════════════════════

void FastKMV::update(const char* seq, uint64_t length) {
    const int K = kmer_size_;
    if (length < static_cast<uint64_t>(K)) return;

    const uint64_t loc_seed = seed_;

    // ── Phase 0: SIMD bulk sequence encoding ────────────────────────────
    std::unique_ptr<uint8_t[]> enc_storage(new uint8_t[length]);
    uint8_t* enc = enc_storage.get();

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
        enc[p] = fkmv_enc_or_invalid((uint8_t)seq[p]);

    // ── Phase 1: ntHash rolling hash ────────────────────────────────────
    uint64_t sf_k[4], sc_ror1[4], sc_km1[4];
    const unsigned km1 = static_cast<unsigned>(K > 1 ? K - 1 : 63);
    for (int i = 0; i < 4; ++i) {
        sf_k[i]    = rol64(NT_SEED_FWD[i], static_cast<unsigned>(K));
        sc_ror1[i] = rol64(NT_SEED_RC[i],  63u);
        sc_km1[i]  = rol64(NT_SEED_RC[i],  km1);
    }

    uint64_t fwd_h = 0, rc_h = 0;
    int inv = 0;
    for (int k = 0; k < K; ++k) {
        uint8_t e = enc[k];
        if (e > 3) { ++inv; e = 0; }
        fwd_h ^= rol64(NT_SEED_FWD[e], static_cast<unsigned>(K - 1 - k));
        rc_h  ^= rol64(NT_SEED_RC[e],  static_cast<unsigned>(k));
    }

    const int      lanes   = 8;
    const uint64_t N_body  = length - static_cast<uint64_t>(K);
    const uint64_t N_batch = (N_body >= 1) ? (N_body / lanes) * lanes : 0;

#if defined(__AVX512F__) && defined(__AVX512DQ__)
    const __m512i vs  = _mm512_set1_epi64((int64_t)loc_seed);
    const __m512i vc1 = _mm512_set1_epi64(0xff51afd7ed558ccdLL);
    const __m512i vc2 = _mm512_set1_epi64(0xc4ceb9fe1a85ec53LL);
    __m512i vthresh = _mm512_set1_epi64((int64_t)threshold_);
#elif defined(__AVX2__)
    const __m256i vs  = _mm256_set1_epi64x((long long)loc_seed);
    const __m256i vc1 = _mm256_set1_epi64x(0xff51afd7ed558ccdLL);
    const __m256i vc2 = _mm256_set1_epi64x(0xc4ceb9fe1a85ec53LL);
    const __m256i vsign = _mm256_set1_epi64x((long long)0x8000000000000000ULL);
    __m256i vthresh = _mm256_set1_epi64x((long long)threshold_);
#endif

    for (uint64_t i = 0; i < N_batch; i += lanes) {
        uint64_t resv[8];
        bool     lane_valid[8];

        for (int j = 0; j < lanes; ++j) {
            const uint64_t pos = i + j;

            lane_valid[j] = (inv == 0);
            resv[j] = lane_valid[j] ? (fwd_h < rc_h ? fwd_h : rc_h) : 0;

            const uint8_t e_out = enc[pos];
            const uint8_t e_in  = enc[pos + K];
            if (e_out > 3) --inv;
            if (e_in  > 3) ++inv;
            const uint8_t oi = (e_out <= 3) ? e_out : 0;
            const uint8_t ii = (e_in  <= 3) ? e_in  : 0;

            fwd_h = rol64(fwd_h, 1)  ^ sf_k[oi]   ^ NT_SEED_FWD[ii];
            rc_h  = rol64(rc_h, 63u) ^ sc_ror1[oi] ^ sc_km1[ii];
        }

#if defined(__AVX512F__) && defined(__AVX512DQ__)
        __mmask8 valid_mask = 0;
        for (int j = 0; j < lanes; ++j)
            if (lane_valid[j]) valid_mask |= (1u << j);
        if (!valid_mask) continue;

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

#if FASTKMV_FMUX_ROUNDS >= 2
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
#elif defined(__AVX2__)
        for (int half = 0; half < 2; ++half) {
            const int base = half * 4;
            if (!(lane_valid[base]|lane_valid[base+1]|lane_valid[base+2]|lane_valid[base+3]))
                continue;

            __m256i vb = _mm256_loadu_si256((const __m256i*)(resv + base));
            __m256i va = _mm256_xor_si256(vb, vs);
            __m256i vt = _mm256_srli_epi64(va, 33);
            vb = _mm256_xor_si256(va, vt);
            va = avx2_mullo_epi64(vb, vc1);
            vt = _mm256_srli_epi64(va, 33);
            vb = _mm256_xor_si256(va, vt);
            va = avx2_mullo_epi64(vb, vc2);
            vt = _mm256_srli_epi64(va, 33);
            vb = _mm256_xor_si256(va, vt);

#ifndef PMH_FAST_HASH
            va = _mm256_xor_si256(vb, vs);
            vt = _mm256_srli_epi64(va, 33);
            vb = _mm256_xor_si256(va, vt);
            va = avx2_mullo_epi64(vb, vc1);
            vt = _mm256_srli_epi64(va, 33);
            vb = _mm256_xor_si256(va, vt);
            va = avx2_mullo_epi64(vb, vc2);
            vt = _mm256_srli_epi64(va, 33);
            vb = _mm256_xor_si256(va, vt);
#endif

            __m256i vkeys = _mm256_srli_epi64(vb, 11);
            __m256i cmp = _mm256_cmpgt_epi64(
                _mm256_xor_si256(vthresh, vsign),
                _mm256_xor_si256(vkeys,   vsign));
            int pass_bits = _mm256_movemask_pd(_mm256_castsi256_pd(cmp));

            if (pass_bits) {
                uint64_t keyv[4];
                _mm256_storeu_si256((__m256i*)keyv, vkeys);
                for (int j = 0; j < 4; ++j) {
                    if ((pass_bits & (1 << j)) && lane_valid[base + j])
                        insertKey(keyv[j]);
                }
                vthresh = _mm256_set1_epi64x((long long)threshold_);
            }
        }
#else
        if (!(lane_valid[0]|lane_valid[1]|lane_valid[2]|lane_valid[3]|
              lane_valid[4]|lane_valid[5]|lane_valid[6]|lane_valid[7])) continue;

        for (int j = 0; j < lanes; ++j) {
            if (!lane_valid[j]) continue;
            uint64_t h0 = mc::murmur3_fmix(resv[j], loc_seed);
#if FASTKMV_FMUX_ROUNDS >= 2
            uint64_t h1 = mc::murmur3_fmix(h0, loc_seed);
#else
            uint64_t h1 = h0;
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
            uint64_t canon = (fwd_h < rc_h) ? fwd_h : rc_h;
            uint64_t h0 = mc::murmur3_fmix(canon, loc_seed);
#if FASTKMV_FMUX_ROUNDS >= 2
            uint64_t h1 = mc::murmur3_fmix(h0, loc_seed);
#else
            uint64_t h1 = h0;
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
            rc_h  = rol64(rc_h, 63u) ^ sc_ror1[oi] ^ sc_km1[ii];
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// jaccard  –  standard KMV bottom-k intersection estimator
// ═══════════════════════════════════════════════════════════════════════════

double FastKMV::jaccard(const FastKMV& other) const {
    assert(k_ == other.k_);
    ensureSorted();
    other.ensureSorted();

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

FastKMV FastKMV::merge(const FastKMV& other) const {
    assert(k_ == other.k_ && kmer_size_ == other.kmer_size_);
    ensureSorted();
    other.ensureSorted();

    FastKMV ret(k_, kmer_size_, seed_);

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
    ret.sorted_ = true;
    return ret;
}

// ═══════════════════════════════════════════════════════════════════════════
// printSketch
// ═══════════════════════════════════════════════════════════════════════════

void FastKMV::printSketch() const {
    ensureSorted();
    std::fprintf(stdout,
                 "FastKMV k=%u kmer=%d fill=%u/%u vals[0..19]: ",
                 k_, kmer_size_, size_, k_);
    for (uint32_t i = 0; i < k_ && i < 20; ++i)
        std::fprintf(stdout, "%016lx ", (unsigned long)vals_[i]);
    if (k_ > 20) std::fprintf(stdout, "...");
    std::fprintf(stdout, "\n");
}

#undef FKMV_VALID
