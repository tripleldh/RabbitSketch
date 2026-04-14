/**
 * SetSketch extreme-optimized implementation.
 *
 * add hot path: 1 integer compare with cached global_thresh_ → skip 95%+ hashes.
 * update():     SIMD hash + inline add with global + per-register early exit.
 * cardinality:  SIMD gather from precomputed base_inv_pow_ table.
 * union_size:   inline max + SIMD gather, zero allocation.
 */
#include "SetSketch.h"
#include "Sketch.h"
#include "MurmurHash3.h"
#include "hash_int.h"
#include <immintrin.h>
#include <algorithm>
#include <cstring>
#include <cstdio>

using namespace Sketch;

// ── SIMD equal-register counter ───────────────────────────────────────────────
namespace {
static int setsketch_count_equal_regs(const uint8_t* __restrict__ a,
                                       const uint8_t* __restrict__ b,
                                       int n) {
  int count = 0, i = 0;
#if defined(__AVX512BW__)
  for (; i + 64 <= n; i += 64) {
    __m512i va = _mm512_loadu_si512((const void*)(a + i));
    __m512i vb = _mm512_loadu_si512((const void*)(b + i));
    count += (int)__builtin_popcountll(
        (uint64_t)_mm512_cmpeq_epi8_mask(va, vb));
  }
#endif
#if defined(__AVX2__)
  for (; i + 32 <= n; i += 32) {
    __m256i va = _mm256_loadu_si256((const __m256i*)(a + i));
    __m256i vb = _mm256_loadu_si256((const __m256i*)(b + i));
    __m256i eq = _mm256_cmpeq_epi8(va, vb);
    count += (int)__builtin_popcount((uint32_t)_mm256_movemask_epi8(eq));
  }
#endif
  for (; i < n; i++)
    count += (int)(a[i] == b[i]);
  return count;
}

static inline double setsketch_sum_max_registers(const uint8_t* __restrict__ c1,
                                                 const uint8_t* __restrict__ c2,
                                                 int m,
                                                 const double* __restrict__ baseInvPow) {
  double sum = 0.0;
  int i = 0;
#if defined(__AVX512BW__) && defined(__AVX512F__)
  __m512d vsum0 = _mm512_setzero_pd();
  __m512d vsum1 = _mm512_setzero_pd();
  for (; i + 16 <= m; i += 16) {
    __m128i va = _mm_loadu_si128((const __m128i*)(c1 + i));
    __m128i vb = _mm_loadu_si128((const __m128i*)(c2 + i));
    __m128i vmax = _mm_max_epu8(va, vb);
    __m256i vidx0 = _mm256_cvtepu8_epi32(vmax);
    __m128i vmax_hi = _mm_srli_si128(vmax, 8);
    __m256i vidx1 = _mm256_cvtepu8_epi32(vmax_hi);
    vsum0 = _mm512_add_pd(vsum0, _mm512_i32gather_pd(vidx0, baseInvPow, 8));
    vsum1 = _mm512_add_pd(vsum1, _mm512_i32gather_pd(vidx1, baseInvPow, 8));
  }
  sum += _mm512_reduce_add_pd(_mm512_add_pd(vsum0, vsum1));
#elif defined(__AVX2__)
  __m256d vacc0 = _mm256_setzero_pd();
  __m256d vacc1 = _mm256_setzero_pd();
  __m256d vacc2 = _mm256_setzero_pd();
  __m256d vacc3 = _mm256_setzero_pd();
  for (; i + 16 <= m; i += 16) {
    __m128i va = _mm_loadu_si128((const __m128i*)(c1 + i));
    __m128i vb = _mm_loadu_si128((const __m128i*)(c2 + i));
    __m128i vmax = _mm_max_epu8(va, vb);

    __m128i idx0_8 = vmax;
    __m128i idx1_8 = _mm_srli_si128(vmax, 8);
    __m256i idx0_32 = _mm256_cvtepu8_epi32(idx0_8);
    __m256i idx1_32 = _mm256_cvtepu8_epi32(idx1_8);
    __m128i idx0_lo = _mm256_castsi256_si128(idx0_32);
    __m128i idx0_hi = _mm256_extracti128_si256(idx0_32, 1);
    __m128i idx1_lo = _mm256_castsi256_si128(idx1_32);
    __m128i idx1_hi = _mm256_extracti128_si256(idx1_32, 1);

    vacc0 = _mm256_add_pd(vacc0, _mm256_i32gather_pd(baseInvPow, idx0_lo, 8));
    vacc1 = _mm256_add_pd(vacc1, _mm256_i32gather_pd(baseInvPow, idx0_hi, 8));
    vacc2 = _mm256_add_pd(vacc2, _mm256_i32gather_pd(baseInvPow, idx1_lo, 8));
    vacc3 = _mm256_add_pd(vacc3, _mm256_i32gather_pd(baseInvPow, idx1_hi, 8));
  }
  __m256d vacc01 = _mm256_add_pd(vacc0, vacc1);
  __m256d vacc23 = _mm256_add_pd(vacc2, vacc3);
  __m256d vacc = _mm256_add_pd(vacc01, vacc23);
  alignas(32) double buf[4];
  _mm256_store_pd(buf, vacc);
  sum += buf[0] + buf[1] + buf[2] + buf[3];
#endif
  for (; i < m; i++) {
    uint8_t r = (c1[i] > c2[i]) ? c1[i] : c2[i];
    sum += baseInvPow[r];
  }
  return sum;
}
} // namespace

// ── Constructor: precompute threshold + base_inv_pow tables ───────────────────
SetSketch::SetSketch(int np, double base, double a)
    : np_(np), q_(62), base_(base), a_(a),
      min_reg_(0), value_(0.0), is_calculated_(0) {
  assert(np >= 4 && np <= 16);
  assert(base > 1.0);
  assert(a > 0.0);

  const uint64_t m = 1ULL << np;
  core_.resize(m, 0);
  shift_ = 64 - np;
  mask_u_ = (shift_ >= 64) ? ~0ULL : (1ULL << shift_) - 1;

  const double log_base = std::log(base);
  factor_ = m * (base - 1.0) / (base * log_base * a);

  // base^(-k) table
  for (int k = 0; k < 64; k++)
    base_inv_pow_[k] = std::pow(base, -(double)k);

  // CDF integer thresholds: threshold[k] = exp(-a * base^(-k)) * 2^shift_
  const double scale = (double)(mask_u_ + 1);
  for (int k = 0; k <= (int)q_; k++) {
    double cdf = std::exp(-a * base_inv_pow_[k]);
    uint64_t t = (uint64_t)(cdf * scale);
    thresholds_[k] = (t > mask_u_) ? mask_u_ : t;
  }
  thresholds_[q_ + 1] = mask_u_ + 1; // sentinel

  // Global lower bound
  count_at_min_ = (uint32_t)m;
  global_thresh_ = thresholds_[0];
}

// ── Recompute min_reg_ by scanning all registers ─────────────────────────────
void SetSketch::recompute_min() {
  uint8_t mn = core_[0];
  for (size_t i = 1; i < core_.size(); i++)
    if (core_[i] < mn) mn = core_[i];
  min_reg_ = mn;
  uint32_t cnt = 0;
  for (size_t i = 0; i < core_.size(); i++)
    cnt += (core_[i] == mn);
  count_at_min_ = cnt;
  global_thresh_ = thresholds_[mn];
}

// ── add_slow: fallback used by remainder loop (not hot path) ──────────────────
void SetSketch::add_slow(uint64_t hashval) {
  const uint64_t rest = hashval & mask_u_;
  if (rest < global_thresh_) return;

  const uint32_t index = (uint32_t)(hashval >> shift_);
  const uint8_t cur = core_[index];
  if (rest < thresholds_[cur]) return;

  uint8_t k = cur + 1;
  while (k < q_ && rest >= thresholds_[k]) ++k;
  core_[index] = k;
  is_calculated_ = 0;

  if (cur == min_reg_) {
    if (--count_at_min_ == 0)
      recompute_min();
  }
}

// ── Cardinality: SIMD gather from base_inv_pow_ table ─────────────────────────
void SetSketch::ensure_cardinality() const {
  if (is_calculated_) return;
  const size_t sz = core_.size();
  const double* __restrict__ tbl = base_inv_pow_;
  const uint8_t* __restrict__ c = core_.data();
  double sum = 0.0;

#if defined(__AVX512F__)
  size_t i = 0;
  __m512d vsum0 = _mm512_setzero_pd();
  __m512d vsum1 = _mm512_setzero_pd();
  for (; i + 16 <= sz; i += 16) {
    // 8 elements per gather
    __m128i vidx8_0 = _mm_loadl_epi64((__m128i*)(c + i));
    __m256i vidx32_0 = _mm256_cvtepu8_epi32(vidx8_0);
    __m512d vals0 = _mm512_i32gather_pd(vidx32_0, tbl, 8);
    vsum0 = _mm512_add_pd(vsum0, vals0);

    __m128i vidx8_1 = _mm_loadl_epi64((__m128i*)(c + i + 8));
    __m256i vidx32_1 = _mm256_cvtepu8_epi32(vidx8_1);
    __m512d vals1 = _mm512_i32gather_pd(vidx32_1, tbl, 8);
    vsum1 = _mm512_add_pd(vsum1, vals1);
  }
  sum = _mm512_reduce_add_pd(_mm512_add_pd(vsum0, vsum1));
  for (; i < sz; i++) sum += tbl[c[i]];
#else
  for (size_t i = 0; i < sz; i++) sum += tbl[c[i]];
#endif

  value_ = (sum > 1e-300) ? factor_ / sum : 0.0;
  is_calculated_ = 1;
}

double SetSketch::cardinality() const {
  ensure_cardinality();
  return value_;
}

// ── union_size: inline max + SIMD gather, no allocation ───────────────────────
double SetSketch::union_size(const SetSketch& other) const {
  const int sz = static_cast<int>(core_.size());
  const double* __restrict__ tbl = base_inv_pow_;
  const uint8_t* __restrict__ c1 = core_.data();
  const uint8_t* __restrict__ c2 = other.core_.data();
  const double sum = setsketch_sum_max_registers(c1, c2, sz, tbl);

  return (sum > 1e-300) ? factor_ / sum : 0.0;
}

double SetSketch::jaccard_index(const SetSketch& other) const {
  double us = union_size(other);
  if (us <= 0.0) return 0.0;
  double c1 = cardinality();
  double c2 = other.cardinality();
  double inter = c1 + c2 - us;
  return (inter > 0.0) ? inter / us : 0.0;
}

SetSketch SetSketch::merge(const SetSketch& other) const {
  assert(np_ == other.np_ && base_ == other.base_ && a_ == other.a_);
  SetSketch ret(np_, base_, a_);
  for (size_t i = 0; i < core_.size(); ++i)
    ret.core_[i] = std::max(core_[i], other.core_[i]);
  ret.is_calculated_ = 0;
  ret.recompute_min();
  return ret;
}

double SetSketch::equalRegisterFraction(const SetSketch& other) const {
  const int n = static_cast<int>(core_.size());
  if (n == 0 || n != static_cast<int>(other.core_.size())) return 0.0;
  return (double)setsketch_count_equal_regs(core_.data(), other.core_.data(), n) / n;
}

double SetSketch::distanceFiltered(const SetSketch& other,
                                   double min_jaccard,
                                   double prefilter_factor) const {
  if (equalRegisterFraction(other) < min_jaccard * prefilter_factor)
    return -1.0;
  return distance(other);
}

void SetSketch::printSketch() {
  fprintf(stdout, "SetSketch core[%zu]: ", core_.size());
  for (size_t i = 0; i < core_.size() && i < 20; i++)
    fprintf(stdout, "%u ", core_[i]);
  if (core_.size() > 20) fprintf(stdout, "...");
  fprintf(stdout, "\n");
}

// ── update(seq): rolling k-mer + SIMD hash + INLINE add with global filter ───
static const uint8_t ENCODE_LUT[256] = {
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
#define ENC(c)   (ENCODE_LUT[(uint8_t)(c)])
#define COMP(e)  ((uint8_t)((e) <= 3 ? 3 - (e) : 255))
#define VALID(e) ((e) <= 3)

void SetSketch::update(char* seq) {
  update(seq, strlen(seq));
}

void SetSketch::update(char* seq, size_t len) {
  const uint64_t LENGTH = static_cast<uint64_t>(len);
  const int KMERLEN = 32;
  if (LENGTH < (uint64_t)KMERLEN) return;

  // Cache class members in locals for the hot loop
  const uint32_t shift = shift_;
  const uint64_t mask = mask_u_;
  const uint32_t qmax = q_;
  uint8_t*       core = core_.data();
  const uint64_t* thresh = thresholds_;
  uint8_t  loc_min_reg     = min_reg_;
  uint32_t loc_count_at_min = count_at_min_;
  uint64_t loc_global_thresh = global_thresh_;
  uint8_t  loc_is_calc = is_calculated_;
  bool     loc_min_dirty = false;

  uint64_t fwd_enc = 0, rev_enc = 0;
  int invalid_count = 0;
  for (int k = 0; k < KMERLEN; k++) {
    uint8_t ef = ENC(seq[k]);
    if (!VALID(ef)) invalid_count++;
    fwd_enc = (fwd_enc << 2) | (VALID(ef) ? (ef & 3u) : 0u);
    uint8_t er = VALID(ef) ? (COMP(ef) & 3u) : 0u;
    rev_enc = (rev_enc >> 2) | (static_cast<uint64_t>(er) << (2 * (KMERLEN - 1)));
  }

  const int lanes = 8;
  const uint64_t total_kmers = LENGTH - KMERLEN + 1;
  const uint64_t N = (total_kmers / lanes) * lanes;

  for (uint64_t i = 0; i < N; i += lanes) {
    uint64_t resv[8];
    bool lane_valid[8];
    for (int j = 0; j < lanes; j++) {
      lane_valid[j] = (invalid_count == 0);
      resv[j] = lane_valid[j] ? ((fwd_enc <= rev_enc) ? fwd_enc : rev_enc) : 0;

      uint8_t ef_out = ENC(seq[i + j]);
      uint8_t ef_in  = ENC(seq[i + j + KMERLEN]);
      if (!VALID(ef_out)) invalid_count--;
      if (!VALID(ef_in))  invalid_count++;
      fwd_enc = (fwd_enc << 2) | (VALID(ef_in) ? (ef_in & 3u) : 0u);
      uint8_t er_in = VALID(ef_in) ? (COMP(ef_in) & 3u) : 0u;
      rev_enc = (rev_enc >> 2) | ((uint64_t)er_in << (2 * (KMERLEN - 1)));
    }

    // ── SIMD hash (identical to HLL) ──────────────────────────────────────
    uint64_t hashvalv[8];
#if defined(__AVX512F__) && defined(__AVX512DQ__)
    __m512i vb = _mm512_loadu_si512((void*)resv);
    __m512i vseed = _mm512_set1_epi64(42);
    __m512i va = _mm512_xor_epi64(vb, vseed);
    __m512i vtmp = _mm512_srli_epi64(va, 33);
    vb = _mm512_xor_epi64(va, vtmp);
    va = _mm512_mullo_epi64(vb, _mm512_set1_epi64(0xff51afd7ed558ccdULL));
    vtmp = _mm512_srli_epi64(va, 33);
    vb = _mm512_xor_epi64(va, vtmp);
    va = _mm512_mullo_epi64(vb, _mm512_set1_epi64(0xc4ceb9fe1a85ec53ULL));
    vtmp = _mm512_srli_epi64(va, 33);
    vb = _mm512_xor_epi64(va, vtmp);
    _mm512_storeu_si512(hashvalv, vb);
#else
    for (int j = 0; j < lanes; j++)
      hashvalv[j] = mc::murmur3_fmix(resv[j], 42);
#endif

    // ── SIMD index/rest extraction + global filter ────────────────────────
#if defined(__AVX512F__)
    {
      __m512i vhash = _mm512_loadu_si512((void*)hashvalv);
      __m512i vmask = _mm512_set1_epi64(mask);
      __m512i vrest = _mm512_and_epi64(vhash, vmask);

      // Global early exit: compare all 8 rests against global threshold
      __m512i vgt = _mm512_set1_epi64(loc_global_thresh);
      // pass_mask bit j=1 → rest[j] >= global_thresh → might update
      __mmask8 pass_mask = _mm512_cmpge_epu64_mask(vrest, vgt);

      // Quick check: if all 8 are filtered, skip entirely
      if (pass_mask != 0) {
        uint64_t restv[8];
        _mm512_storeu_si512(restv, vrest);

        for (int j = 0; j < lanes; j++) {
          if (!lane_valid[j]) continue;
          if (!(pass_mask & (1 << j))) continue;  // global filter

          uint32_t idx = (uint32_t)(hashvalv[j] >> shift);
          uint64_t rest = restv[j];
          uint8_t cur = core[idx];

          // Per-register early exit
          if (rest < thresh[cur]) continue;

          // Find new register value (scan thresholds, ~1-2 steps)
          uint8_t k = cur + 1;
          while (k < qmax && rest >= thresh[k]) ++k;

          core[idx] = k;
          loc_is_calc = 0;

          // Update global lower-bound tracking
          if (cur == loc_min_reg) {
            if (!loc_min_dirty && loc_count_at_min > 0 && --loc_count_at_min == 0)
              loc_min_dirty = true;
          }
        }
      }
    }
#else
    // Scalar path with global filter
    for (int j = 0; j < lanes; j++) {
      if (!lane_valid[j]) continue;
      uint64_t rest = hashvalv[j] & mask;
      if (rest < loc_global_thresh) continue;  // global filter

      uint32_t idx = (uint32_t)(hashvalv[j] >> shift);
      uint8_t cur = core[idx];
      if (rest < thresh[cur]) continue;  // per-register filter

      uint8_t k = cur + 1;
      while (k < qmax && rest >= thresh[k]) ++k;
      core[idx] = k;
      loc_is_calc = 0;

      if (cur == loc_min_reg) {
        if (!loc_min_dirty && loc_count_at_min > 0 && --loc_count_at_min == 0)
          loc_min_dirty = true;
      }
    }
#endif
  }

  // ── Remainder loop ──────────────────────────────────────────────────────
  for (uint64_t i = N; i < total_kmers; ++i) {
    if (invalid_count == 0) {
      uint64_t res = (fwd_enc <= rev_enc) ? fwd_enc : rev_enc;
      uint64_t hashval = mc::murmur3_fmix(res, 42);
      uint64_t rest = hashval & mask;
      if (rest >= loc_global_thresh) {
        uint32_t idx = (uint32_t)(hashval >> shift);
        uint8_t cur = core[idx];
        if (rest >= thresh[cur]) {
          uint8_t k = cur + 1;
          while (k < qmax && rest >= thresh[k]) ++k;
          core[idx] = k;
          loc_is_calc = 0;
          if (cur == loc_min_reg) {
            if (!loc_min_dirty && loc_count_at_min > 0 && --loc_count_at_min == 0)
              loc_min_dirty = true;
          }
        }
      }
    }
    uint8_t ef_out = ENC(seq[i]);
    uint8_t ef_in  = ENC(seq[i + KMERLEN]);
    if (!VALID(ef_out)) invalid_count--;
    if (!VALID(ef_in))  invalid_count++;
    fwd_enc = (fwd_enc << 2) | (VALID(ef_in) ? (ef_in & 3u) : 0u);
    uint8_t er_in = VALID(ef_in) ? (COMP(ef_in) & 3u) : 0u;
    rev_enc = (rev_enc >> 2) | ((uint64_t)er_in << (2 * (KMERLEN - 1)));
  }

  // Write back locals. If min tracking was exhausted, rebuild once at end.
  if (loc_min_dirty) {
    is_calculated_ = loc_is_calc;
    recompute_min();
  } else {
    min_reg_ = loc_min_reg;
    count_at_min_ = loc_count_at_min;
    global_thresh_ = loc_global_thresh;
    is_calculated_ = loc_is_calc;
  }
}

// ── containment ──────────────────────────────────────────────────────────────
// C(this ⊆ other) = (|A| + |B| - |AUB|) / |A|
double SetSketch::containment(const SetSketch& other) const
{
  const double card_a = cardinality();
  if (card_a <= 0.0) return 0.0;
  const double card_b = other.cardinality();
  const double us     = union_size(other);
  const double inter  = card_a + card_b - us;
  return (inter > 0.0) ? inter / card_a : 0.0;
}

// ── ani ──────────────────────────────────────────────────────────────────────
// ANI = (2J / (1+J))^(1/kmer_size)   (Mash / Ondov et al. 2016)
double SetSketch::ani(const SetSketch& other, int kmer_size) const
{
  const double j = jaccard_index(other);
  if (j <= 0.0) return 0.0;
  if (j >= 1.0) return 1.0;
  return std::pow(2.0 * j / (1.0 + j), 1.0 / static_cast<double>(kmer_size));
}

// ── inverted index: block key extraction ─────────────────────────────────────
void SetSketch::getBlockKeys(std::vector<uint32_t>& keys, int blockSize) const {
    const int m = getM();
    const int nBlocks = m / blockSize;
    keys.resize(nBlocks);
    for (int b = 0; b < nBlocks; b++) {
        keys[b] = blockHash(static_cast<uint32_t>(b),
                            core_[b * 3],
                            core_[b * 3 + 1],
                            core_[b * 3 + 2]);
    }
}

// ── inverted index: exact Jaccard from flat core arrays (SIMD) ──────────────
double SetSketch::jaccardFromCores(
    const uint8_t* __restrict__ c1,
    const uint8_t* __restrict__ c2,
    int m,
    const double* __restrict__ baseInvPow,
    double factor,
    double card1, double card2)
{
    const double sum = setsketch_sum_max_registers(c1, c2, m, baseInvPow);
    if (sum <= 1e-300) return 0.0;
    double us    = factor / sum;
    double inter = card1 + card2 - us;
    return (inter > 0.0) ? inter / us : 0.0;
}

double SetSketch::jaccardFromCoresEarlyAbort(
    const uint8_t* __restrict__ c1,
    const uint8_t* __restrict__ c2,
    int m,
    const double* __restrict__ baseInvPow,
    double factor,
    double card1, double card2,
    double minJaccard)
{
    if (m <= 0) return 0.0;
    if (minJaccard <= 0.0) {
        return jaccardFromCores(c1, c2, m, baseInvPow, factor, card1, card2);
    }
    const double cardSum = card1 + card2;
    const double maxTerm = baseInvPow[0];
    double sum = 0.0;
    for (int i = 0; i < m; ++i) {
        uint8_t r = (c1[i] > c2[i]) ? c1[i] : c2[i];
        sum += baseInvPow[r];

        const int remain = m - i - 1;
        const double maxPossibleSum = sum + maxTerm * static_cast<double>(remain);
        const double maxPossibleJ = (cardSum * maxPossibleSum / factor) - 1.0;
        if (maxPossibleJ < minJaccard) return -1.0;
    }
    if (sum <= 1e-300) return 0.0;
    const double us = factor / sum;
    const double inter = cardSum - us;
    return (inter > 0.0) ? inter / us : 0.0;
}

#undef ENC
#undef COMP
#undef VALID
