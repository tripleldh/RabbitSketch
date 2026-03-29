/**
 * ProbMinHash4 sketch for probability Jaccard similarity estimation.
 *
 * Based on: Otmar Ertl, "ProbMinHash – A Class of Locality-Sensitive Hash
 * Algorithms for the (Probability) Jaccard Similarity", IEEE TKDE 2020.
 * https://arxiv.org/abs/1911.00675
 *
 * Optimizations vs. baseline:
 *  1. Ziggurat exponential distribution  (~5x faster than -log(u))
 *  2. SIMD batched murmur3_fmix hashing  (AVX-512: 8-lane; AVX2: 8-lane emulated)
 *  3. Early return in addHashFromRng BEFORE perm_.reset()
 *  4. Cached class members in hot loop locals
 *  5. Global-max pre-filter to skip most addHashFromRng calls
 *  6. SIMD Jaccard / merge  (AVX-512 / AVX2 / scalar)
 *  7. PermStream version overflow protection
 *  8. No strlen(): caller passes length directly (eliminates one full scan)
 *  9. TED parameters packed into TedParam[] (better cache locality)
 * 10. tracker_ leaves ARE the sketch registers (no duplicate regs_ storage)
 * 11. O(m) build_from_leaves() for copy / merge (vs. O(m log m) on-line updates)
 * 12. cur_max loaded once per batch (register); refreshed only after real update
 * 13. bool[8] lane_valid; all-N batch skipped via branchless OR check
 * 14. Optional single-fmix fast path: compile with -DPMH_FAST_HASH
 *
 * Interface mirrors HyperLogLog / SetSketch in this library:
 *   update(seq, length) – ingest a DNA/RNA sequence (k-mer rolling hash)
 *   distance(other)     – probability Jaccard distance  [0, 1]
 *   jaccard(other)      – probability Jaccard similarity [0, 1]
 *   merge(other)        – element-wise min merge (returns new sketch)
 */

#ifndef _PROBMH_H_
#define _PROBMH_H_

#include <cstdint>
#include <cmath>
#include <limits>
#include <memory>
#include <cassert>

namespace Sketch {

// ── Internal: tournament-tree max-value tracker ───────────────────────────
// Maintains the current maximum over m registers; isUpdatePossible() lets
// ProbMinHash4 skip elements that can no longer improve any register.
// Leaf values (indices 0..m-1) double as the sketch registers themselves,
// eliminating the need for a separate regs_ array.
class ProbMHMaxTracker {
public:
    explicit ProbMHMaxTracker(uint32_t m);
    void   reset(double infinity);
    bool   update(uint32_t idx, double value);
    bool   isUpdatePossible(double value) const;
    double getMax() const;

    // Direct access to the m leaf values (= sketch registers).
    double*       leaves()       noexcept { return v_.get(); }
    const double* leaves() const noexcept { return v_.get(); }

    // O(m) bottom-up rebuild of internal nodes from already-set leaves.
    // Call after bulk-writing leaves (copy / merge) instead of m individual
    // update() calls, which would cost O(m log m).
    void build_from_leaves();

    friend void swap(ProbMHMaxTracker& a, ProbMHMaxTracker& b) noexcept;

private:
    uint32_t m_;
    uint32_t lastIdx_;
    std::unique_ptr<double[]> v_;   // [0..m-1] leaves, [m..2m-2] internal nodes
};

void swap(ProbMHMaxTracker& a, ProbMHMaxTracker& b) noexcept;

// ── Internal: lazy Fisher-Yates permutation stream ────────────────────────
// SoA (val_ + ver_arr_) for better cache locality when reading version in next().
class ProbMHPermStream {
public:
    explicit ProbMHPermStream(uint32_t m);
    void     reset();
    uint32_t next(uint64_t& rng_state);   // advances rng_state in-place
    friend void swap(ProbMHPermStream& a, ProbMHPermStream& b) noexcept;

private:
    uint32_t m_;
    uint32_t idx_;
    uint32_t ver_;
    std::unique_ptr<uint32_t[]> val_;     // logical value at position
    std::unique_ptr<uint32_t[]> ver_arr_; // version stamp
};

void swap(ProbMHPermStream& a, ProbMHPermStream& b) noexcept;

// ── ProbMinHash4 sketch ────────────────────────────────────────────────────
class ProbMinHash4 {
public:
    /**
     * @param m         number of sketch registers (sketch size); must be > 1
     * @param kmer_size k-mer length for sequence hashing (default 21)
     * @param seed      64-bit seed for wyhash-based RNG (default 42)
     * @param max_L     max register updates per element (Route C truncation).
     *                  0 = unlimited (original ProbMinHash4 behaviour).
     *                  1/2/4/… = Top-L truncation (faster, slight accuracy loss).
     */
    explicit ProbMinHash4(uint32_t m = 1024,
                          int      kmer_size = 21,
                          uint64_t seed = 42,
                          uint32_t max_L = 0);

    ~ProbMinHash4() = default;
    ProbMinHash4(const ProbMinHash4&);
    ProbMinHash4& operator=(ProbMinHash4 other);
    ProbMinHash4(ProbMinHash4&&) = default;
    ProbMinHash4& operator=(ProbMinHash4&&) = default;

    /**
     * Ingest a DNA/RNA sequence.
     * Canonical k-mers (forward vs reverse-complement minimum) are used.
     * Can be called multiple times to add more sequences to the same sketch.
     * The caller must supply the sequence length (avoids an extra strlen pass).
     */
    void update(const char* seq, uint64_t length);

    /**
     * Return the probability Jaccard similarity in [0, 1].
     */
    double jaccard(const ProbMinHash4& other) const;

    /**
     * Return probability Jaccard distance = 1 - jaccard().
     */
    double distance(const ProbMinHash4& other) const {
        return 1.0 - jaccard(other);
    }

    /**
     * Merge two sketches by taking element-wise minimum of hash values.
     * Both sketches must have the same parameters.
     */
    ProbMinHash4 merge(const ProbMinHash4& other) const;

    /**
     * Return pointer to the m raw register values (for serialization / LSH).
     * The pointer remains valid as long as the sketch is alive and unmodified.
     */
    const double* getRegisters() const noexcept { return tracker_.leaves(); }

    int      getKmerSize()  const { return kmer_size_; }
    uint32_t getM()         const { return m_; }
    uint32_t getMaxL()      const { return max_L_; }

    void printSketch() const;

private:
    void addHash(uint64_t h);
    void addHashFromRng(uint64_t rng);

    // Packed TED (Truncated-Exponential Distribution) parameters.
    // Replaces five separate arrays (boundaries_, ted_rate_, ted_c1/c2/c3_)
    // for better spatial locality: hot-loop accesses stride over one array.
    struct TedParam {
        double boundary;  // normalised cumulative boundary; [0] = 1.0
        double gap;       // boundary[i] - boundary[i-1]; used for delta (i >= 1)
        double c1, c2, c3;
    };

    uint32_t           m_;
    int                kmer_size_;
    uint64_t           seed_;
    uint32_t           max_L_;    // Route C: max updates per element (m_ = unlimited)

    std::unique_ptr<TedParam[]> ted_params_;     // [m-1]
    double                      firstBoundaryInv_;

    ProbMHMaxTracker   tracker_;
    ProbMHPermStream   perm_;
};

// ── One-Permutation ProbMinHash (Route A) ──────────────────────────────────
//
// Each element maps to exactly ONE bucket via b(x) = hash1(x) % m, with a
// single key g(x) = Uniform(0,1) from hash2(x).  Each register keeps the
// minimum key:  S_j = min_{x: b(x)=j} g(x).
//
// Compared to ProbMinHash4:
//   • No PermStream (eliminated perm_.reset() / perm_.next())
//   • No MaxTracker tournament tree (flat array + scalar global_max)
//   • No TED sampling; one Uniform(0,1) per element
//   • O(1) per-element update (vs O(log m) tree + permutation)
//   • global_max pre-filter rejects ~99%+ elements after warm-up
//
// Trade-offs:
//   • Empty buckets possible when #elements < m  (rare for genomic data)
//   • Jaccard estimator uses count(A[j]==B[j] && finite) / m
//   • Statistical distribution differs from original ProbMinHash
//
class ProbMinHash4OP {
public:
    explicit ProbMinHash4OP(uint32_t m = 1024,
                            int      kmer_size = 21,
                            uint64_t seed = 42);

    ~ProbMinHash4OP() = default;
    ProbMinHash4OP(const ProbMinHash4OP&);
    ProbMinHash4OP& operator=(ProbMinHash4OP other);
    ProbMinHash4OP(ProbMinHash4OP&&) = default;
    ProbMinHash4OP& operator=(ProbMinHash4OP&&) = default;

    void update(const char* seq, uint64_t length);

    double jaccard(const ProbMinHash4OP& other) const;

    double distance(const ProbMinHash4OP& other) const {
        return 1.0 - jaccard(other);
    }

    ProbMinHash4OP merge(const ProbMinHash4OP& other) const;

    const double* getRegisters() const noexcept { return regs_.get(); }
    uint32_t getM()        const { return m_; }
    int      getKmerSize() const { return kmer_size_; }
    uint32_t numEmpty()    const;

    void printSketch() const;

private:
    void addHash(uint64_t h);
    void refreshMax();

    uint32_t m_;
    int      kmer_size_;
    uint64_t seed_;

    std::unique_ptr<double[]> regs_;
    double   global_max_;
    uint32_t num_nonempty_;
};

// ── ProbKMV (Route B) ──────────────────────────────────────────────────────
//
// For each element x, generate a single key g(x) = Uniform(0,1) from
// hash(x), then keep the k smallest keys across all elements:
//
//     sketch(A) = bottom-k { g(x) : x ∈ A }
//
// This is a KMV (K Minimum Values) structure: canonical k-mers are hashed with
// MurmurHash3_x64_128 (8-byte key, lower 64 bits → [0,1); same family as MinHash).
// Optional -DPMH_FAST_HASH uses one round instead of two.  Jaccard uses the
// standard KMV two-pointer merge:
//
//     J ≈ |S_A ∩ S_B in bottom-k of S_A ∪ S_B| / k
//
// Compared to ProbMinHash4 / ProbMinHash4OP:
//   • Structure is a sorted uint64_t[k] array (bottom-k integer keys)
//   • Keys are 53-bit unsigned integers (h >> 11), order-equivalent to
//     the [0,1) double representation; all comparisons stay in integer
//     domain (no FP conversions in insert, jaccard, or threshold filter)
//   • No per-bucket partitioning; no empty-bucket problem
//   • Merge = sorted merge of two bottom-k lists → take k smallest
//   • Jaccard = intersection-over-union in merged bottom-k
//   • One hash per element (no bucket hash needed)
//   • Naturally deduplicates repeated k-mers
//
// Trade-offs:
//   • Warmup: O(1) append into 2k buffer; first compactify sorts,
//     deduplicates, and truncates to sorted bottom-k.
//   • Steady state: sorted array with O(log k) binary_search + O(k)
//     memmove — identical to original, hardware-optimised path.
//   • Threshold pre-filter rejects >99% of elements after warm-up
//   • This is a KMV sketch, not ProbMinHash proper
//
class ProbKMV {
public:
    static constexpr uint64_t KEY_MAX = (UINT64_C(1) << 53) - 1;

    explicit ProbKMV(uint32_t k = 1024,
                     int      kmer_size = 21,
                     uint64_t seed = 42);

    ~ProbKMV() = default;
    ProbKMV(const ProbKMV&);
    ProbKMV& operator=(ProbKMV other);
    ProbKMV(ProbKMV&&) = default;
    ProbKMV& operator=(ProbKMV&&) = default;

    void update(const char* seq, uint64_t length);

    /**
     * KMV Jaccard estimator: merge the two sorted bottom-k lists, take
     * the k smallest distinct values, and count how many are shared.
     */
    double jaccard(const ProbKMV& other) const;

    double distance(const ProbKMV& other) const {
        return 1.0 - jaccard(other);
    }

    ProbKMV merge(const ProbKMV& other) const;

    const uint64_t* getRegisters() const { ensureSorted(); return vals_.get(); }
    uint32_t getK()        const { return k_; }
    uint32_t getM()        const { return k_; }
    int      getKmerSize() const { return kmer_size_; }
    uint32_t size()        const { ensureSorted(); return size_; }

    void printSketch() const;

private:
    void addHash(uint64_t h);
    void insertKey(uint64_t key);
    void compactify() const;
    void ensureSorted() const;

    uint32_t k_;
    int      kmer_size_;
    uint64_t seed_;
    uint32_t buf_cap_;

    std::unique_ptr<uint64_t[]> vals_;
    mutable uint32_t size_;
    mutable uint64_t threshold_;
    mutable bool     sorted_;
};

} // namespace Sketch

#endif // _PROBMH_H_
