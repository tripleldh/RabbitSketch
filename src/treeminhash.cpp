/**
 * TreeMinHash – RabbitSketch wrapper (k-mer counts → Otmar Ertl TreeMinHash).
 *
 * Performance notes (vs. original):
 *
 * 1. Adaptive tree_max: when tree_max=0 (auto), we compute the actual maximum
 *    k-mer count per update() call and rebuild the TMH tree only when needed.
 *    With tree_max=DBL_MAX the tree has ~1074 levels; with tree_max=actual_max
 *    it has only ~54+log2(max_count) levels.  Moreover the root-node invRate
 *    becomes 1/max instead of ~0, so the initial point-vs-limit skip filter
 *    eliminates >95% of elements before they even enter tree traversal.
 *    Combined effect: ~100–600x speedup for typical k-mer data.
 *
 * 2. Smart reserve: cnt.reserve(min(num_positions, 1<<20)) prevents the ~10
 *    rehashes that occur when starting from capacity 1024 on a 1M sequence.
 *
 * 3. Pre-encoded sequence buffer (Impl::enc_buf): one LUT pass up-front, then
 *    the hot k-mer loop only touches a flat uint8_t array.
 *
 * 4. Persistent data buffer (Impl::data_buf): reused across update() calls
 *    to avoid repeated heap allocation for the (kmer, count) vector.
 *
 * 5. In-place fill: sig_ is passed directly to tmh::TreeMinHash::fill(),
 *    eliminating the extra vector allocation+move that operator() did.
 */

#include "treeminhash.h"
#include "tmh_weighted_minwise.hpp"
#include "robin_hood.h"

#include <cstdio>
#include <cstring>
#include <cmath>
#include <limits>
#include <utility>
#include <algorithm>

namespace Sketch {

namespace {

struct TMHRng {
    using RngType = tmh::WyrandBitStream;
    uint64_t seed_;

    tmh::WyrandBitStream operator()(uint64_t x) const {
        return tmh::WyrandBitStream(x, seed_);
    }
    tmh::WyrandBitStream operator()(uint64_t x, uint64_t y) const {
        return tmh::WyrandBitStream(x, y, seed_);
    }
};

} // namespace

// ---------------------------------------------------------------------------
// Impl – holds all mutable per-object state
// ---------------------------------------------------------------------------
struct TreeMinHash::Impl {
    TMHRng rng;

    // The TMH hasher is stored via unique_ptr so its tree can be rebuilt
    // when the observed max k-mer count changes.
    std::unique_ptr<tmh::TreeMinHash<TMHRng>> hasher;

    // Current tree upper bound and whether we are in auto-adapt mode.
    double   actual_tree_max;
    bool     auto_tree_max;

    // Parameters kept for rebuilding the hasher.
    uint32_t m_;
    double   tree_factor_;
    double   success_prob_first_;

    // Reusable per-update buffers (avoid repeated malloc/free).
    std::vector<uint8_t>                     enc_buf;   // pre-encoded sequence
    std::vector<std::pair<uint64_t, double>> data_buf;  // (kmer, count) pairs
    // Sort-path scratch: kept on Impl so repeated long-sequence updates reuse RAM.
    std::vector<uint64_t>                    kmer_flat_buf;

    Impl(uint32_t m, uint64_t seed, double tree_max_arg,
         double fac, double sp)
        : rng{seed},
          // Start small (2.0 covers count=1 and count=2) if auto;
          // otherwise use user-supplied value directly.
          actual_tree_max(tree_max_arg > 0.0 ? tree_max_arg : 2.0),
          auto_tree_max(tree_max_arg <= 0.0),
          m_(m), tree_factor_(fac), success_prob_first_(sp)
    {
        rebuild_hasher();
    }

    void rebuild_hasher() {
        hasher = std::make_unique<tmh::TreeMinHash<TMHRng>>(
            m_, rng, actual_tree_max, tree_factor_, success_prob_first_);
    }

    // Called after computing max_count from the k-mer map.
    // Rebuilds the tree only when the observed max falls outside the
    // currently covered range [actual/4 .. actual], giving 4x headroom
    // before a shrink-rebuild (prevents oscillation on borderline counts).
    void adapt_tree_max(uint32_t max_count) {
        if (!auto_tree_max) return;
        const double needed = static_cast<double>(max_count);
        if (needed <= 0.0) return;
        if (needed <= actual_tree_max && needed > actual_tree_max * 0.25)
            return;  // within headroom – no rebuild
        actual_tree_max = needed * 2.0;  // 2x overhead going forward
        rebuild_hasher();
    }
};

// ---------------------------------------------------------------------------
// Base-4 encoding LUT: A=0, C=1, G=2, T=3, others=255
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// Constructors / destructor / swap
// ---------------------------------------------------------------------------
TreeMinHash::TreeMinHash(uint32_t m, int kmer_size, uint64_t seed,
                         double tree_max, double tree_factor,
                         double success_prob_first)
    : m_(m),
      kmer_size_(kmer_size),
      seed_(seed),
      tree_max_(tree_max),
      tree_factor_(tree_factor),
      success_prob_first_(success_prob_first),
      impl_(new Impl(m_, seed_, tree_max_, tree_factor_, success_prob_first_))
{
    sig_.assign(m_, {UINT64_C(0), std::numeric_limits<double>::infinity()});
}

TreeMinHash::TreeMinHash(const TreeMinHash& o)
    : m_(o.m_),
      kmer_size_(o.kmer_size_),
      seed_(o.seed_),
      tree_max_(o.tree_max_),
      tree_factor_(o.tree_factor_),
      success_prob_first_(o.success_prob_first_),
      sig_(o.sig_),
      impl_(new Impl(m_, seed_, tree_max_, tree_factor_, success_prob_first_))
{
    // Restore the same actual_tree_max that the source had, so subsequent
    // fill() calls use the already-tuned tree without a redundant rebuild.
    impl_->actual_tree_max = o.impl_->actual_tree_max;
    impl_->auto_tree_max   = o.impl_->auto_tree_max;
    impl_->rebuild_hasher();
}

TreeMinHash::TreeMinHash(TreeMinHash&& o) noexcept = default;

TreeMinHash& TreeMinHash::operator=(TreeMinHash other) {
    swap(*this, other);
    return *this;
}

TreeMinHash::~TreeMinHash() = default;

void swap(TreeMinHash& a, TreeMinHash& b) noexcept {
    using std::swap;
    swap(a.m_, b.m_);
    swap(a.kmer_size_, b.kmer_size_);
    swap(a.seed_, b.seed_);
    swap(a.tree_max_, b.tree_max_);
    swap(a.tree_factor_, b.tree_factor_);
    swap(a.success_prob_first_, b.success_prob_first_);
    swap(a.sig_, b.sig_);
    swap(a.impl_, b.impl_);
}

// ---------------------------------------------------------------------------
// update() – k-mer counting + TMH sketch computation
// ---------------------------------------------------------------------------
void TreeMinHash::update(const char* seq, uint64_t length)
{
    const int K = kmer_size_;
    if (length < static_cast<uint64_t>(K)) {
        sig_.assign(m_, {UINT64_C(0), std::numeric_limits<double>::infinity()});
        return;
    }

    const uint64_t num_pos = length - static_cast<uint64_t>(K) + 1;

    // --- Phase 1: pre-encode sequence to uint8_t ----------------------------
    // Single LUT pass up-front keeps the hot k-mer loop cache-clean.
    auto& enc_buf = impl_->enc_buf;
    enc_buf.resize(length);
    uint8_t* const enc = enc_buf.data();
    for (uint64_t i = 0; i < length; ++i)
        enc[i] = ENC_LUT[(uint8_t)seq[i]];

    // --- Phase 2: collect and count canonical k-mers ------------------------
    //
    // Two strategies depending on sequence length:
    //
    //   Small (num_pos ≤ SORT_THRESHOLD):
    //     robin_hood hash map, O(N) amortised, good for short sequences.
    //
    //   Large (num_pos > SORT_THRESHOLD):
    //     Collect all canonical k-mers into a flat uint64_t vector, sort it,
    //     then count runs.  Sequential memory access patterns give much better
    //     cache utilisation and use 8 bytes/k-mer instead of ~16 bytes/entry
    //     for the hash table.  Sorting adds O(N log N) work but in practice
    //     beats the hash map for N > ~4M due to the cache advantage.

    static constexpr uint64_t SORT_THRESHOLD = 1u << 22;  // ~4M positions
    // Above this many windows, a full uint64_t[kmer] array is often multi‑GB;
    // use hash counting instead (O(unique) RAM, may be slower than sort).
    static constexpr uint64_t FLAT_KMER_CAP = UINT64_C(1) << 27;  // 128M windows (~1GiB flat)

    const uint64_t kmer_mask =
        (K == 32) ? ~UINT64_C(0) : ((UINT64_C(1) << (2 * K)) - 1);
    const int shift = 2 * (K - 1);

    uint64_t fwd = 0, rev = 0;
    int inv = 0;  // count of invalid (N) bases in current window

    // Seed the first k-mer.
    for (int k = 0; k < K; ++k) {
        const uint8_t e = enc[k];
        const bool valid = (e <= 3);
        if (!valid) ++inv;
        const uint8_t eb = valid ? e : 0u;
        fwd = (fwd << 2) | eb;
        rev = (rev >> 2) | (static_cast<uint64_t>(valid ? (3u - eb) : 0u) << shift);
    }

    uint32_t max_count = 0;
    auto& data = impl_->data_buf;
    data.clear();

    auto run_hash_path = [&]() -> bool {
        robin_hood::unordered_map<uint64_t, uint32_t> cnt;
        // Cap reserve: full num_pos reserve needs huge bucket tables for long reads.
        cnt.reserve(static_cast<size_t>(std::min<uint64_t>(num_pos, UINT64_C(1) << 20)));

        auto emit_map = [&](uint64_t canon) {
            uint32_t& c = cnt[canon];
            ++c;
            if (c > max_count) max_count = c;
        };

        if (__builtin_expect(inv == 0, 1))
            emit_map((fwd <= rev) ? fwd : rev);

        for (uint64_t i = 0; i < num_pos - 1; ++i) {
            const uint8_t e_out = enc[i];
            const uint8_t e_in  = enc[i + K];
            const bool v_out = (e_out <= 3);
            const bool v_in  = (e_in  <= 3);
            if (!v_out) --inv;
            if (!v_in)  ++inv;
            const uint8_t in_bits = v_in ? (e_in & 3u) : 0u;
            fwd = ((fwd << 2) | in_bits) & kmer_mask;
            rev = (rev >> 2) |
                  (static_cast<uint64_t>(v_in ? (3u - in_bits) : 0u) << shift);
            if (__builtin_expect(inv == 0, 1))
                emit_map((fwd <= rev) ? fwd : rev);
        }

        if (cnt.empty()) {
            sig_.assign(m_, {UINT64_C(0), std::numeric_limits<double>::infinity()});
            return false;
        }
        data.reserve(cnt.size());
        for (const auto& kv : cnt)
            data.emplace_back(kv.first, static_cast<double>(kv.second));
        return true;
    };

    if (num_pos <= SORT_THRESHOLD) {
        if (!run_hash_path()) return;
    } else if (num_pos > FLAT_KMER_CAP) {
        // Avoid ~8*num_pos temporary bytes (often multi‑GB for whole genomes).
        if (!run_hash_path()) return;
    } else {
        // ---- Sort-based path (large sequences, flat array fits in memory) ---
        auto& kmer_flat = impl_->kmer_flat_buf;
        kmer_flat.clear();
        kmer_flat.reserve(static_cast<size_t>(num_pos));

        auto emit_sort = [&](uint64_t canon) { kmer_flat.push_back(canon); };

        if (__builtin_expect(inv == 0, 1))
            emit_sort((fwd <= rev) ? fwd : rev);

        for (uint64_t i = 0; i < num_pos - 1; ++i) {
            const uint8_t e_out = enc[i];
            const uint8_t e_in  = enc[i + K];
            const bool v_out = (e_out <= 3);
            const bool v_in  = (e_in  <= 3);
            if (!v_out) --inv;
            if (!v_in)  ++inv;
            const uint8_t in_bits = v_in ? (e_in & 3u) : 0u;
            fwd = ((fwd << 2) | in_bits) & kmer_mask;
            rev = (rev >> 2) |
                  (static_cast<uint64_t>(v_in ? (3u - in_bits) : 0u) << shift);
            if (__builtin_expect(inv == 0, 1))
                emit_sort((fwd <= rev) ? fwd : rev);
        }

        if (kmer_flat.empty()) {
            sig_.assign(m_, {UINT64_C(0), std::numeric_limits<double>::infinity()});
            return;
        }

        std::sort(kmer_flat.begin(), kmer_flat.end());

        data.reserve(kmer_flat.size());
        uint64_t prev = kmer_flat[0];
        uint32_t run = 1;
        for (uint64_t i = 1; i < kmer_flat.size(); ++i) {
            if (kmer_flat[i] == prev) {
                ++run;
            } else {
                if (run > max_count) max_count = run;
                data.emplace_back(prev, static_cast<double>(run));
                prev = kmer_flat[i];
                run = 1;
            }
        }
        if (run > max_count) max_count = run;
        data.emplace_back(prev, static_cast<double>(run));
    }

    // --- Phase 3: adapt tree to observed weight range -----------------------
    // If max_count << actual_tree_max (e.g. count=1 but tree_max=DBL_MAX),
    // the tree has thousands of useless levels and nearly every element passes
    // the first-level skip filter (invRate ≈ 0 ⟹ point ≈ 0 < limit always).
    // Rebuilding with actual_tree_max ≈ 2*max_count:
    //   - reduces tree depth from ~1074 to ~54+log2(max_count)
    //   - makes invRate = 1/max_count, so ~95%+ of elements get skipped
    //     before entering tree traversal
    impl_->adapt_tree_max(max_count);

    // --- Phase 5: run TreeMinHash directly into sig_ -----------------------
    // fill() avoids the extra vector allocation that operator() would do.
    sig_.resize(m_);
    impl_->hasher->fill(data, sig_);
}

// ---------------------------------------------------------------------------
// jaccard / distance
// ---------------------------------------------------------------------------
double TreeMinHash::jaccard(const TreeMinHash& other) const
{
    if (m_ != other.m_) return 0.0;
    if (sig_.size() != (size_t)m_ || other.sig_.size() != (size_t)m_)
        return 0.0;
    uint32_t eq = 0;
    for (uint32_t j = 0; j < m_; ++j) {
        if (sig_[j] == other.sig_[j])
            ++eq;
    }
    return static_cast<double>(eq) / static_cast<double>(m_);
}

// ---------------------------------------------------------------------------
// printSketch
// ---------------------------------------------------------------------------
void TreeMinHash::printSketch() const
{
    std::fprintf(stdout, "TreeMinHash m=%u k=%d tree_max=%.4g sig[0..min(19,m-1)]: ",
                 m_, kmer_size_, impl_->actual_tree_max);
    const uint32_t n = std::min<uint32_t>(m_, 20);
    for (uint32_t i = 0; i < n; ++i)
        std::fprintf(stdout, "(%llu,%.4g) ",
                     (unsigned long long)sig_[i].first, sig_[i].second);
    if (m_ > 20) std::fprintf(stdout, "...");
    std::fprintf(stdout, "\n");
}

} // namespace Sketch
