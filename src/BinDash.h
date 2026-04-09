/**
 * BinDash – b-bit One Permutation Hashing sketch for Jaccard similarity.
 *
 * Algorithm based on: Zhao, X. "BinDash, software for fast genome distance
 * estimation on a typical personal laptop", Bioinformatics 2019.
 *
 * Key property: extremely low memory per sketch.  Default parameters
 * (sketchsize64=32, bbits=16) yield 2048 bins packed into 512 uint64_t
 * words = 4096 bytes per sketch, compared to 16384 bytes for a traditional
 * 2048-register MinHash storing full 64-bit hashes (4x compression).
 *
 * Interface mirrors HyperLogLog / ProbMinHash4 / FastKMV in this library:
 *   update(seq, length) – ingest a DNA/RNA sequence (canonical k-mers)
 *   jaccard(other)      – estimated Jaccard similarity [0, 1]
 *   distance(other)     – 1 - jaccard
 */

#ifndef _BINDASH_H_
#define _BINDASH_H_

#include <cstdint>
#include <vector>
#include <cstring>

namespace Sketch {

class BinDash {
public:
    /**
     * @param sketchsize64  number of 64-bin groups (total bins = sketchsize64 * 64)
     * @param kmer_size     k-mer length for sequence hashing (default 21)
     * @param bbits         bits retained per bin minimum (default 16)
     * @param seed          hash seed (default 42)
     */
    explicit BinDash(uint32_t sketchsize64 = 32,
                     int      kmer_size    = 21,
                     uint32_t bbits        = 16,
                     uint64_t seed         = 42);

    ~BinDash() = default;
    BinDash(const BinDash&) = default;
    BinDash& operator=(const BinDash&) = default;
    BinDash(BinDash&&) = default;
    BinDash& operator=(BinDash&&) = default;

    /**
     * Ingest a DNA/RNA sequence.  Can be called multiple times to add more
     * sequences (e.g. multi-contig genomes).  The sketch is automatically
     * finalized (densified + bit-packed) on first jaccard()/distance() call.
     */
    void update(const char* seq, uint64_t length);
    void update(const char* seq);

    /**
     * Explicitly finalize (densify + bit-pack + free temporary signs_).
     * Called automatically on first jaccard()/distance(), but should be
     * called explicitly before multi-threaded distance computation.
     */
    void finalize();

    /**
     * Estimated Jaccard similarity in [0, 1].
     * Triggers finalization on first call.
     */
    double jaccard(const BinDash& other) const;

    double distance(const BinDash& other) const {
        return 1.0 - jaccard(other);
    }

    const uint64_t* getSignatures() const { return usigs_.data(); }
    uint32_t getNumWords()  const { return sketchsize64_ * bbits_; }
    uint32_t getNumBins()   const { return nbins_; }
    uint32_t getBbits()     const { return bbits_; }
    int      getKmerSize()  const { return kmer_size_; }
    size_t   memoryBytes()  const { return usigs_.size() * sizeof(uint64_t); }

    /**
     * Extract the truncated bbits-bit value for a single bin (inverse of packBits).
     * Must be called after finalize().
     */
    uint16_t getBinValue(uint32_t bin) const {
        const uint32_t bb   = bbits_;
        const uint32_t word = bin / 64;
        const uint32_t bit  = bin % 64;
        uint16_t val = 0;
        for (uint32_t b = 0; b < bb; b++)
            val |= static_cast<uint16_t>(((usigs_[word * bb + b] >> bit) & 1ULL) << b);
        return val;
    }

private:
    uint32_t sketchsize64_;
    int      kmer_size_;
    uint32_t bbits_;
    uint64_t seed_;
    uint32_t nbins_;

    std::vector<uint64_t> signs_;   // [nbins_] per-bin minimum hash (temporary)
    std::vector<uint64_t> usigs_;   // [sketchsize64_*bbits_] bit-packed (final)
    mutable bool finalized_;

    void ensureFinalized() const;
    void densify();
    void packBits();
};

} // namespace Sketch

#endif // _BINDASH_H_
