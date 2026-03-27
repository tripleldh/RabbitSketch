/**
 * TreeMinHash – weighted Jaccard sketch (Otmar Ertl, treeminhash).
 *
 * Builds a multiset of k-mers from DNA/RNA (canonical k-mer, counts as weights),
 * then runs the TreeMinHash algorithm on (kmer_hash, count) pairs.
 *
 * @see https://github.com/oertl/treeminhash
 */

#ifndef _RABBITSKETCH_TREEMINHASH_H_
#define _RABBITSKETCH_TREEMINHASH_H_

#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

namespace Sketch {

class TreeMinHash {
public:
    /**
     * @param m                  signature length (sketch size)
     * @param kmer_size          k-mer length (1–32)
     * @param seed               RNG seed for WyrandBitStream
     * @param tree_max           upper bound on k-mer weights; 0 = auto-adapt
     *                           per update() call (recommended, avoids 1000+
     *                           tree levels when weights are small integers)
     * @param tree_factor        geometric split factor in (0,1), default 0.5
     * @param success_prob_first first-pass success probability (default 0.9)
     */
    explicit TreeMinHash(uint32_t m = 1024,
                         int      kmer_size = 21,
                         uint64_t seed = 42,
                         double   tree_max = 0.0,
                         double   tree_factor = 0.5,
                         double   success_prob_first = 0.9);

    TreeMinHash(const TreeMinHash&);
    TreeMinHash& operator=(TreeMinHash other);
    TreeMinHash(TreeMinHash&&) noexcept;
    TreeMinHash& operator=(TreeMinHash&&) noexcept;

    ~TreeMinHash();

    void update(const char* seq, uint64_t length);

    double jaccard(const TreeMinHash& other) const;

    double distance(const TreeMinHash& other) const { return 1.0 - jaccard(other); }

    uint32_t getM() const { return m_; }
    int      getKmerSize() const { return kmer_size_; }
    uint64_t getSeed() const { return seed_; }

    const std::vector<std::pair<uint64_t, double>>& getSignature() const { return sig_; }

    void printSketch() const;

private:
    uint32_t m_;
    int      kmer_size_;
    uint64_t seed_;
    double   tree_max_;
    double   tree_factor_;
    double   success_prob_first_;

    std::vector<std::pair<uint64_t, double>> sig_;

    struct Impl;
    std::unique_ptr<Impl> impl_;

    friend void swap(TreeMinHash& a, TreeMinHash& b) noexcept;
};

} // namespace Sketch

#endif
