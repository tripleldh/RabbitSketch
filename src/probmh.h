/**
 * ProbMinHash4 sketch for probability Jaccard similarity estimation.
 *
 * Based on: Otmar Ertl, "ProbMinHash – A Class of Locality-Sensitive Hash
 * Algorithms for the (Probability) Jaccard Similarity", IEEE TKDE 2020.
 * https://arxiv.org/abs/1911.00675
 *
 * Interface mirrors HyperLogLog / SetSketch in this library:
 *   update(char* seq)  – ingest a DNA/RNA sequence (k-mer rolling hash)
 *   distance(other)    – probability Jaccard distance  [0, 1]
 *   jaccard(other)     – probability Jaccard similarity [0, 1]
 *   merge(other)       – element-wise min merge (returns new sketch)
 */

#ifndef _PROBMH_H_
#define _PROBMH_H_

#include <cstdint>
#include <vector>
#include <cmath>
#include <limits>
#include <memory>
#include <cassert>

namespace Sketch {

// ── Internal: tournament-tree max-value tracker ───────────────────────────
// Maintains the current maximum over m registers; isUpdatePossible() lets
// ProbMinHash4 skip elements that can no longer improve any register.
class ProbMHMaxTracker {
public:
    explicit ProbMHMaxTracker(uint32_t m);
    void   reset(double infinity);
    bool   update(uint32_t idx, double value);
    bool   isUpdatePossible(double value) const;
    double getMax() const;
    friend void swap(ProbMHMaxTracker& a, ProbMHMaxTracker& b) noexcept;

private:
    uint32_t m_;
    uint32_t lastIdx_;
    std::unique_ptr<double[]> v_;
};

void swap(ProbMHMaxTracker& a, ProbMHMaxTracker& b) noexcept;

// ── Internal: lazy Fisher-Yates permutation stream ────────────────────────
// Produces the next element of a random permutation of {0,..,m-1} on demand.
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
    std::unique_ptr<std::pair<uint32_t,uint32_t>[]> pv_;
};

void swap(ProbMHPermStream& a, ProbMHPermStream& b) noexcept;

// ── ProbMinHash4 sketch ────────────────────────────────────────────────────
class ProbMinHash4 {
public:
    /**
     * @param m         number of sketch registers (sketch size); must be > 1
     * @param kmer_size k-mer length for sequence hashing (default 21)
     * @param seed      64-bit seed for wyhash-based RNG (default 42)
     */
    explicit ProbMinHash4(uint32_t m = 1024,
                          int      kmer_size = 21,
                          uint64_t seed = 42);

    ~ProbMinHash4() = default;
    ProbMinHash4(const ProbMinHash4&);
    ProbMinHash4& operator=(ProbMinHash4 other);
    ProbMinHash4(ProbMinHash4&&) = default;
    ProbMinHash4& operator=(ProbMinHash4&&) = default;

    /**
     * Ingest a DNA/RNA sequence.
     * Canonical k-mers (forward vs reverse-complement minimum) are used.
     * Can be called multiple times to add more sequences to the same sketch.
     */
    void update(char* seq);

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
     * Return the raw register values (for serialization / external use).
     */
    const std::vector<double>& getRegisters() const { return regs_; }

    int      getKmerSize()  const { return kmer_size_; }
    uint32_t getM()         const { return m_; }

    void printSketch() const;

private:
    // Feed a single 64-bit canonical k-mer hash into ProbMinHash4 update.
    void addHash(uint64_t h);

    uint32_t           m_;
    int                kmer_size_;
    uint64_t           seed_;

    std::vector<double> regs_;          // m_ hash registers (all inf initially)

    // ProbMinHash4 precomputed tables
    std::unique_ptr<double[]>  boundaries_;          // [m-1]  normalised boundaries
    std::unique_ptr<double[]>  ted_rate_;            // [m-1]  truncated-exp rates
    std::unique_ptr<double[]>  ted_c1_;              // [m-1]  (exp(r)-1)/r
    std::unique_ptr<double[]>  ted_c2_;              // [m-1]
    std::unique_ptr<double[]>  ted_c3_;              // [m-1]
    double                     firstBoundaryInv_;    // 1 / log1p(1/(m-1))

    ProbMHMaxTracker   tracker_;
    ProbMHPermStream   perm_;
};

} // namespace Sketch

#endif // _PROBMH_H_
