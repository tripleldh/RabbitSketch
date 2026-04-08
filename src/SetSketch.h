/**
 * SetSketch - Extreme-optimized: O(1) integer-only add with global lower-bound filter,
 * SIMD-accelerated update(), and SIMD-gather cardinality/distance.
 *
 * Hot path: 1 comparison against cached global threshold → skip >95% of hashes.
 * Register distribution: P(K ≤ k) = exp(-a * base^(-k)).
 */
#ifndef _SETSKETCH_H_
#define _SETSKETCH_H_

#include <cstdint>
#include <vector>
#include <cmath>
#include <cassert>

namespace Sketch {

class SetSketch {
public:
  SetSketch(int np, double base = 2.0, double a = 20.0);
  ~SetSketch() = default;

  void update(char* seq);
  SetSketch merge(const SetSketch& other) const;
  double cardinality() const;
  double jaccard_index(const SetSketch& other) const;
  double distance(const SetSketch& other) const { return 1.0 - jaccard_index(other); }

  /// Containment of *this in other: |A ∩ B| / |A|.
  ///   C(A⊆B) = (|A| + |B| - |AUB|) / |A|
  double containment(const SetSketch& other) const;

  /// Average Nucleotide Identity from Jaccard similarity.
  ///   ANI = (2J / (1+J))^(1/kmer_size)
  /// @param kmer_size  k-mer length used during sketching (default 32).
  double ani(const SetSketch& other, int kmer_size = 32) const;
  const std::vector<uint8_t>& getCore() const { return core_; }
  double equalRegisterFraction(const SetSketch& other) const;
  double distanceFiltered(const SetSketch& other,
                          double min_jaccard,
                          double prefilter_factor = 0.3) const;
  void printSketch();

  // Expose internals for inline distance computation in test harness
  double getFactor() const { return factor_; }
  const double* getBaseInvPow() const { return base_inv_pow_; }
  int getM() const { return (int)(1ULL << np_); }

private:
  void add_slow(uint64_t hashval);
  void recompute_min();
  double union_size(const SetSketch& other) const;
  void ensure_cardinality() const;

  std::vector<uint8_t> core_;
  uint32_t np_;
  uint32_t q_;
  double base_;
  double a_;
  double factor_;

  // Precomputed lookup tables
  uint64_t thresholds_[64];            // CDF thresholds (stack-allocated, always in cache)
  double   base_inv_pow_[64];          // base^(-k) for k=0..63
  uint64_t mask_u_;
  uint32_t shift_;

  // Global lower-bound tracking (like HLL's implicit filter via clz)
  uint8_t  min_reg_;                   // current min register value
  uint32_t count_at_min_;              // registers still at min_reg_
  uint64_t global_thresh_;             // = thresholds_[min_reg_], cached in register

  mutable double value_;
  mutable uint8_t is_calculated_;
};

} // namespace Sketch

#endif
