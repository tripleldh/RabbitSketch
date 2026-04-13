/**
 * FastKMV – fast K Minimum Values sketch for genomic k-mers.
 *
 * ntHash rolling hash + murmur3 fmix; SIMD encoding and batched fmix where
 * available.  Compile with -DFASTKMV_NO_FMUX to use raw ntHash>>11 keys (ablation).
 * Not part of ProbMinHash; standalone KMV / bottom-k estimator.
 *
 * Jaccard: standard KMV two-pointer merge on sorted bottom-k keys.
 */

#ifndef _FASTKMV_H_
#define _FASTKMV_H_

#include <cstdint>
#include <memory>
#include <cassert>

namespace Sketch {

class FastKMV {
public:
    static constexpr uint64_t KEY_MAX = (UINT64_C(1) << 53) - 1;

    explicit FastKMV(uint32_t k = 1024,
                     int      kmer_size = 21,
                     uint64_t seed = 42);

    ~FastKMV() = default;
    FastKMV(const FastKMV&);
    FastKMV& operator=(FastKMV other);
    FastKMV(FastKMV&&) = default;
    FastKMV& operator=(FastKMV&&) = default;

    void update(const char* seq, uint64_t length);

    /**
     * KMV Jaccard: merge two sorted bottom-k lists (conceptually), count
     * overlap in the k smallest distinct values.
     */
    double jaccard(const FastKMV& other) const;

    double distance(const FastKMV& other) const {
        return 1.0 - jaccard(other);
    }

    /**
     * KMV cardinality estimate: (k-1) * KEY_MAX / tau_k.
     * Returns the exact count when the sketch is not yet full.
     */
    double cardinality() const;

    /**
     * Containment of *this in other: |A ∩ B| / |A|.
     * Uses the cardinality-based formula:
     *   C(A⊆B) = J * (|A| + |B|) / (|A| * (1 + J))
     */
    double containment(const FastKMV& other) const;

    /**
     * Average Nucleotide Identity estimated from Jaccard similarity.
     * ANI = (2J / (1+J))^(1/k)  where k is the k-mer size stored at build time.
     */
    double ani(const FastKMV& other) const;

    FastKMV merge(const FastKMV& other) const;

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
    mutable std::unique_ptr<uint8_t[]> enc_buf_;
    mutable uint64_t enc_cap_;
    mutable uint32_t size_;
    mutable uint64_t threshold_;
    mutable bool     sorted_;
};

} // namespace Sketch

#endif // _FASTKMV_H_
