/**
 * ProbMinHash4OP – One-Permutation ProbMinHash implementation.
 *
 * Route A: compress the original "per-element permutation stream" into a
 * single bucket assignment b(x) and a single key g(x):
 *
 *     S_j = min_{x : b(x) = j}  g(x)
 *
 * where  b(x) = hash1(x) % m        (uniform bucket)
 *        g(x) = Uniform(0,1)         (from hash2(x))
 *
 * Key design decisions:
 *  1. Flat double[m] registers (no tournament tree)
 *  2. global_max_ pre-filter: after all registers are filled, any element
 *     whose key >= max(regs) is skipped in O(1).  For typical genomic data
 *     (millions of k-mers, m=1024), this rejects >99% of elements.
 *  3. global_max_ is recomputed via O(m) scan only when the register that
 *     WAS the global max gets updated (rare event).
 *  4. SIMD batched murmur3_fmix for k-mer hashing (AVX-512 / scalar fallback)
 *  5. SIMD Jaccard comparison and merge (AVX-512 / AVX2 / scalar)
 *  6. Uniform(0,1) keys instead of Exp(1): for the equality-based Jaccard
 *     estimator, the key distribution does not affect correctness; uniform
 *     avoids the log() / ziggurat cost entirely.
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
// ProbMinHash4OP  –  constructor / copy / assign
// ═══════════════════════════════════════════════════════════════════════════

ProbMinHash4OP::ProbMinHash4OP(uint32_t m, int kmer_size, uint64_t seed)
    : m_(m), kmer_size_(kmer_size), seed_(seed),
      regs_(new double[m]),
      global_max_(std::numeric_limits<double>::infinity()),
      num_nonempty_(0)
{
    assert(m > 1);
    assert(kmer_size >= 1 && kmer_size <= 32);
    std::fill_n(regs_.get(), m, std::numeric_limits<double>::infinity());
}

ProbMinHash4OP::ProbMinHash4OP(const ProbMinHash4OP& o)
    : m_(o.m_), kmer_size_(o.kmer_size_), seed_(o.seed_),
      regs_(new double[o.m_]),
      global_max_(o.global_max_),
      num_nonempty_(o.num_nonempty_)
{
    std::copy(o.regs_.get(), o.regs_.get() + m_, regs_.get());
}

ProbMinHash4OP& ProbMinHash4OP::operator=(ProbMinHash4OP other) {
    std::swap(m_,            other.m_);
    std::swap(kmer_size_,    other.kmer_size_);
    std::swap(seed_,         other.seed_);
    std::swap(regs_,         other.regs_);
    std::swap(global_max_,   other.global_max_);
    std::swap(num_nonempty_, other.num_nonempty_);
    return *this;
}

// ═══════════════════════════════════════════════════════════════════════════
// refreshMax  –  O(m) rescan to find the new global maximum register value
// ═══════════════════════════════════════════════════════════════════════════

void ProbMinHash4OP::refreshMax() {
    double mx = 0;
    const double* r = regs_.get();
    for (uint32_t i = 0; i < m_; ++i)
        if (r[i] > mx) mx = r[i];
    global_max_ = mx;
}

// ═══════════════════════════════════════════════════════════════════════════
// addHash  –  one-permutation insert: bucket + key → array min
//
// Hot path per element:
//   2 × murmur3_fmix  →  1 double conversion  →  1 compare (global_max)
//   →  1 modular bucket  →  1 compare-and-swap on regs_[bucket]
//
// After warm-up, >99% of elements are rejected by the global_max filter.
// ═══════════════════════════════════════════════════════════════════════════

void ProbMinHash4OP::addHash(uint64_t h) {
    uint64_t h1 = mc::murmur3_fmix(h, seed_);
    uint64_t h2 = mc::murmur3_fmix(h1, seed_);

    double key = (double)(h2 >> 11) * PMH_INV2_53;

    if (key >= global_max_) return;

    // Lemire-style fast modular reduction: upper 32 bits of (h1_lo32 * m)
    uint32_t bucket = (uint32_t)(((uint64_t)(uint32_t)h1 * (uint64_t)m_) >> 32);

    double old_val = regs_[bucket];
    if (key < old_val) {
        regs_[bucket] = key;
        const double inf = std::numeric_limits<double>::infinity();
        if (old_val == inf) {
            if (++num_nonempty_ == m_)
                refreshMax();
        } else if (old_val == global_max_) {
            refreshMax();
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// update  –  k-mer rolling hash with SIMD batch hashing
// ═══════════════════════════════════════════════════════════════════════════

void ProbMinHash4OP::update(const char* seq, uint64_t length) {
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

        // ── 8-lane batched murmur3_fmix (AVX-512) for k-mer → h0 ────────
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

        // ── Per-lane: bucket + key + register update ─────────────────────
        double cur_max = global_max_;
        for (int j = 0; j < lanes; ++j) {
            if (!lane_valid[j]) continue;

#ifdef PMH_FAST_HASH
            uint64_t h1 = h0v[j];
#else
            uint64_t h1 = mc::murmur3_fmix(h0v[j], loc_seed);
#endif
            uint64_t h2 = mc::murmur3_fmix(h1, loc_seed);
            double key = (double)(h2 >> 11) * PMH_INV2_53;

            if (key >= cur_max) continue;

            uint32_t bucket = (uint32_t)(((uint64_t)(uint32_t)h1 * (uint64_t)m_) >> 32);
            double old_val = regs_[bucket];

            if (key < old_val) {
                regs_[bucket] = key;
                const double inf = std::numeric_limits<double>::infinity();
                if (old_val == inf) {
                    if (++num_nonempty_ == m_) {
                        refreshMax();
                        cur_max = global_max_;
                    }
                } else if (old_val == cur_max) {
                    refreshMax();
                    cur_max = global_max_;
                }
            }
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
            uint64_t h2 = mc::murmur3_fmix(h1, loc_seed);
            double key = (double)(h2 >> 11) * PMH_INV2_53;

            if (key < global_max_) {
                uint32_t bucket = (uint32_t)(((uint64_t)(uint32_t)h1 * (uint64_t)m_) >> 32);
                double old_val = regs_[bucket];
                if (key < old_val) {
                    regs_[bucket] = key;
                    const double inf = std::numeric_limits<double>::infinity();
                    if (old_val == inf) {
                        if (++num_nonempty_ == m_)
                            refreshMax();
                    } else if (old_val == global_max_) {
                        refreshMax();
                    }
                }
            }
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
// jaccard  (SIMD comparison — AVX-512 / AVX2 / scalar)
// ═══════════════════════════════════════════════════════════════════════════

double ProbMinHash4OP::jaccard(const ProbMinHash4OP& other) const {
    assert(m_ == other.m_);
    const double* __restrict__ a = regs_.get();
    const double* __restrict__ b = other.regs_.get();
    const double inf = std::numeric_limits<double>::infinity();
    int count = 0;
    uint32_t k = 0;

#if defined(__AVX512F__)
    {
        __m512d vinf = _mm512_set1_pd(inf);
        for (; k + 8 <= m_; k += 8) {
            __m512d va = _mm512_loadu_pd(a + k);
            __m512d vb = _mm512_loadu_pd(b + k);
            __mmask8 eq     = _mm512_cmp_pd_mask(va, vb, _CMP_EQ_OQ);
            __mmask8 notinf = _mm512_cmp_pd_mask(va, vinf, _CMP_NEQ_UQ);
            count += __builtin_popcount(eq & notinf);
        }
    }
#elif defined(__AVX2__)
    {
        __m256d vinf = _mm256_set1_pd(inf);
        for (; k + 4 <= m_; k += 4) {
            __m256d va  = _mm256_loadu_pd(a + k);
            __m256d vb  = _mm256_loadu_pd(b + k);
            __m256d ceq = _mm256_cmp_pd(va, vb, _CMP_EQ_OQ);
            __m256d cni = _mm256_cmp_pd(va, vinf, _CMP_NEQ_UQ);
            __m256d res = _mm256_and_pd(ceq, cni);
            count += __builtin_popcount(_mm256_movemask_pd(res));
        }
    }
#endif

    for (; k < m_; ++k)
        count += (a[k] == b[k] && a[k] != inf) ? 1 : 0;

    return static_cast<double>(count) / static_cast<double>(m_);
}

double ProbMinHash4OP::distance(const ProbMinHash4OP& other) const {
    const double j = jaccard(other);
    if (j <= 0.0) return std::numeric_limits<double>::infinity();
    if (j >= 1.0) return 0.0;
    const double ratio = 2.0 * j / (1.0 + j);
    return -std::log(ratio) / static_cast<double>(kmer_size_);
}

// ═══════════════════════════════════════════════════════════════════════════
// merge  –  element-wise min
// ═══════════════════════════════════════════════════════════════════════════

ProbMinHash4OP ProbMinHash4OP::merge(const ProbMinHash4OP& other) const {
    assert(m_ == other.m_ && kmer_size_ == other.kmer_size_);
    ProbMinHash4OP ret(m_, kmer_size_, seed_);

    const double* __restrict__ a = regs_.get();
    const double* __restrict__ b = other.regs_.get();
    double*       __restrict__ r = ret.regs_.get();
    uint32_t k = 0;

#if defined(__AVX512F__)
    for (; k + 8 <= m_; k += 8) {
        __m512d vr = _mm512_min_pd(_mm512_loadu_pd(a + k),
                                   _mm512_loadu_pd(b + k));
        _mm512_storeu_pd(r + k, vr);
    }
#elif defined(__AVX2__)
    for (; k + 4 <= m_; k += 4) {
        __m256d vr = _mm256_min_pd(_mm256_loadu_pd(a + k),
                                   _mm256_loadu_pd(b + k));
        _mm256_storeu_pd(r + k, vr);
    }
#endif
    for (; k < m_; ++k)
        r[k] = std::min(a[k], b[k]);

    const double inf = std::numeric_limits<double>::infinity();
    ret.num_nonempty_ = 0;
    double mx = 0;
    for (uint32_t i = 0; i < m_; ++i) {
        if (r[i] != inf) {
            ret.num_nonempty_++;
            if (r[i] > mx) mx = r[i];
        }
    }
    ret.global_max_ = (ret.num_nonempty_ == m_) ? mx : inf;

    return ret;
}

// ═══════════════════════════════════════════════════════════════════════════
// numEmpty / printSketch
// ═══════════════════════════════════════════════════════════════════════════

uint32_t ProbMinHash4OP::numEmpty() const {
    const double inf = std::numeric_limits<double>::infinity();
    uint32_t cnt = 0;
    for (uint32_t i = 0; i < m_; ++i)
        if (regs_[i] == inf) ++cnt;
    return cnt;
}

void ProbMinHash4OP::printSketch() const {
    std::fprintf(stdout,
                 "ProbMinHash4OP m=%u k=%d empty=%u regs[0..19]: ",
                 m_, kmer_size_, numEmpty());
    for (uint32_t i = 0; i < m_ && i < 20; ++i)
        std::fprintf(stdout, "%.4g ", regs_[i]);
    if (m_ > 20) std::fprintf(stdout, "...");
    std::fprintf(stdout, "\n");
}

#undef PMH_ENC
#undef PMH_COMP
#undef PMH_VALID
