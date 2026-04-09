/**
 * BinDash – b-bit One Permutation Hashing implementation.
 *
 * Algorithm:
 *   1. update(): canonical k-mer sliding window → murmur3_fmix hash →
 *      partition into nbins bins, keep minimum hash per bin.
 *   2. densify(): fill empty bins using optimal densification
 *      (Shrivastava, ICML 2017).
 *   3. packBits(): retain only the lowest bbits of each bin minimum,
 *      bit-pack into sketchsize64 * bbits uint64_t words.
 *   4. jaccard(): XNOR + AND across bbits layers, popcount matching bits.
 */

#include "BinDash.h"
#include "hash_int.h"

#include <algorithm>
#include <cassert>
#include <climits>
#include <cstring>

using namespace Sketch;

// ═══════════════════════════════════════════════════════════════════════════
// DNA encoding (same LUT as ProbMinHash / FastKMV)
// ═══════════════════════════════════════════════════════════════════════════

static const uint8_t BD_ENC_LUT[256] = {
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

#define BD_ENC(c)    (BD_ENC_LUT[(uint8_t)(c)])
#define BD_COMP(e)   ((uint8_t)((e) <= 3 ? 3 - (e) : 255))
#define BD_VALID(e)  ((e) <= 3)

// ═══════════════════════════════════════════════════════════════════════════
// Constructor
// ═══════════════════════════════════════════════════════════════════════════

BinDash::BinDash(uint32_t sketchsize64, int kmer_size,
                 uint32_t bbits, uint64_t seed)
    : sketchsize64_(sketchsize64),
      kmer_size_(kmer_size),
      bbits_(bbits),
      seed_(seed),
      nbins_(sketchsize64 * 64),
      signs_(sketchsize64 * 64, UINT64_MAX),
      usigs_(sketchsize64 * bbits, 0),
      finalized_(false)
{
    assert(sketchsize64 > 0);
    assert(kmer_size >= 1 && kmer_size <= 32);
    assert(bbits >= 1 && bbits <= 64);
}

// ═══════════════════════════════════════════════════════════════════════════
// update – canonical k-mer sliding window + bin-minimum tracking
// ═══════════════════════════════════════════════════════════════════════════

void BinDash::update(const char* seq) {
    update(seq, static_cast<uint64_t>(std::strlen(seq)));
}

void BinDash::update(const char* seq, uint64_t length) {
    finalized_ = false;

    const int K = kmer_size_;
    if (length < static_cast<uint64_t>(K)) return;

    const uint64_t kmer_mask = (K == 32) ? ~0ULL : ((1ULL << (2 * K)) - 1);
    const uint64_t loc_seed  = seed_;
    const uint32_t loc_nbins = nbins_;

    uint64_t fwd = 0, rev = 0;
    int inv = 0;

    for (int k = 0; k < K; ++k) {
        uint8_t ef = BD_ENC(seq[k]);
        if (!BD_VALID(ef)) inv++;
        fwd = (fwd << 2) | (BD_VALID(ef) ? (ef & 3u) : 0u);
        uint8_t er = BD_VALID(ef) ? BD_COMP(ef) : 0u;
        rev = (rev >> 2) | (static_cast<uint64_t>(er & 3u) << (2 * (K - 1)));
    }

    const uint64_t N_body = length - static_cast<uint64_t>(K);

    for (uint64_t i = 0; i <= N_body; ++i) {
        if (inv == 0) {
            const uint64_t canonical = (fwd <= rev) ? fwd : rev;
            uint64_t h = mc::murmur3_fmix(canonical, loc_seed);
            uint32_t binidx = static_cast<uint32_t>((h >> 32) % loc_nbins);
            uint64_t value  = mc::murmur3_fmix(h, 0ULL);
            if (value < signs_[binidx])
                signs_[binidx] = value;
        }

        if (i < N_body) {
            uint8_t ef_out = BD_ENC(seq[i]);
            uint8_t ef_in  = BD_ENC(seq[i + K]);
            if (!BD_VALID(ef_out)) inv--;
            if (!BD_VALID(ef_in))  inv++;
            fwd = ((fwd << 2) | (BD_VALID(ef_in) ? (ef_in & 3u) : 0u)) & kmer_mask;
            uint8_t er_in = BD_VALID(ef_in) ? (BD_COMP(ef_in) & 3u) : 0u;
            rev = (rev >> 2) | (static_cast<uint64_t>(er_in) << (2 * (K - 1)));
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// densify – optimal densification (Shrivastava, ICML 2017)
// ═══════════════════════════════════════════════════════════════════════════

static inline uint64_t univhash2(uint64_t s, uint64_t t) {
    uint64_t x = 1009ULL * s + 1000003ULL * t;
    return (48271ULL * x + 11ULL) % ((1ULL << 31) - 1);
}

void BinDash::densify() {
    uint64_t minval = UINT64_MAX;
    uint64_t maxval = 0;
    for (auto sign : signs_) {
        if (sign < minval) minval = sign;
        if (sign > maxval) maxval = sign;
    }
    if (maxval != UINT64_MAX) return;
    if (minval == UINT64_MAX) return;

    const uint64_t n = signs_.size();
    for (uint64_t i = 0; i < n; i++) {
        uint64_t j = i;
        uint64_t nattempts = 0;
        while (signs_[j] == UINT64_MAX) {
            j = univhash2(i, nattempts) % n;
            nattempts++;
        }
        signs_[i] = signs_[j];
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// packBits – b-bit truncation + column-major bit packing
// ═══════════════════════════════════════════════════════════════════════════

void BinDash::packBits() {
    std::fill(usigs_.begin(), usigs_.end(), 0ULL);
    const uint32_t bb = bbits_;
    for (uint32_t signidx = 0; signidx < nbins_; signidx++) {
        uint64_t sign = signs_[signidx];
        uint32_t leftshift = signidx % 64;
        for (uint32_t b = 0; b < bb; b++) {
            uint64_t bit = (sign >> b) & 1ULL;
            usigs_[signidx / 64 * bb + b] |= (bit << leftshift);
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// ensureFinalized – lazy densify + pack + free temporary signs_
// ═══════════════════════════════════════════════════════════════════════════

void BinDash::finalize() {
    if (finalized_) return;
    densify();
    packBits();
    { std::vector<uint64_t>().swap(signs_); }
    finalized_ = true;
}

void BinDash::ensureFinalized() const {
    if (finalized_) return;
    const_cast<BinDash*>(this)->finalize();
}

// ═══════════════════════════════════════════════════════════════════════════
// jaccard – XNOR + AND across bbits layers, popcount
// ═══════════════════════════════════════════════════════════════════════════

double BinDash::jaccard(const BinDash& other) const {
    assert(sketchsize64_ == other.sketchsize64_);
    assert(bbits_ == other.bbits_);

    ensureFinalized();
    other.ensureFinalized();

    const uint32_t ss64 = sketchsize64_;
    const uint32_t bb   = bbits_;
    uint64_t samebits = 0;

    for (uint32_t i = 0; i < ss64; i++) {
        uint64_t bits = ~0ULL;
        for (uint32_t j = 0; j < bb; j++) {
            bits &= ~(usigs_[i * bb + j] ^ other.usigs_[i * bb + j]);
        }
        samebits += __builtin_popcountll(bits);
    }

    const double maxnbits = static_cast<double>(ss64) * 64.0;
    const double p_match  = static_cast<double>(samebits) / maxnbits;
    const double p_random = 1.0 / static_cast<double>(1ULL << bb);
    double j = (p_match - p_random) / (1.0 - p_random);
    if (j < 0.0) j = 0.0;
    if (j > 1.0) j = 1.0;
    return j;
}

#undef BD_ENC
#undef BD_COMP
#undef BD_VALID
