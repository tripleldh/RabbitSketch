#ifndef _COMMON_H_
#define _COMMON_H_
#include <iostream>
#include <cmath>
#include <cstdint>
#include <vector>
#include <sys/time.h>
#include <unistd.h>
//typedef struct kssd_parameter
//{
//	int half_k;
//	int half_subk;
//	int drlevel;
//	int rev_add_move;
//	int half_outctx_len;
//	int * shuffled_dim;
//	int dim_start;
//	int dim_end;
//	unsigned int kmer_size;
//	int hashSize;
//	int hashLimit;
//	uint64_t domask;
//	uint64_t tupmask;
//	uint64_t undomask0;
//	uint64_t undomask1;
//} kssd_parameter_t;
//
//static const int BaseMap[128] = 
//{
//-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
//-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
//-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
//-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
//-1, 0, -1, 1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, 
//-1, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
//-1, 0, -1, 1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, 
//-1, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
//};
//
static const uint32_t primer[25] = 
{
	251, 509, 1021, 2039, 4093, 8191, 16381,
	32749, 65521, 131071, 262139, 524287,
	1048573, 2097143, 4194301, 8388593, 16777213,
	33554393, 67108859, 134217689, 268435399,
	536870909, 1073741789, 2147483647, 4294967291
};


double get_sec();
unsigned int get_hashSize(int half_k, int drlevel);
uint64_t get_total_system_memory();
int get_progress_bar_size(int total_num);

// ── Per-k-mer Shannon entropy weight (sliding window, O(L)) ──────────────────
//
// w[i] = H(kmer_i) / log2(4)  ∈ [w_min, 1.0]
//
// H = log2(k) - (1/k) * Σ_b n_b * log2(n_b)   (n_b = count of base b in kmer)
//
// Biological intuition:
//   - Low-complexity k-mers (poly-A runs, microsatellites): H ≈ 0 → low weight
//   - High-complexity k-mers (uniform base distribution): H ≈ 2 → high weight
//   - k-mers containing N or invalid bases: weight = w_min
//
// This weight is INTRINSIC to the k-mer sequence (same k-mer → same entropy),
// so it satisfies the ProbMinHash requirement that the weight of a k-mer is
// the same across different sequences.
//
// Parameters:
//   seq   : raw sequence (ACGT or acgt; other chars treated as invalid)
//   L     : sequence length
//   k     : k-mer length
//   w_min : minimum weight (default 0.1; prevents weight=0 for degenerate kmers)
static inline void fill_kmer_entropy_weights(std::vector<double>& w,
                                              const char* seq, uint64_t L, int k,
                                              double w_min = 0.1)
{
    const uint64_t nwin = (L >= (uint64_t)k) ? L - (uint64_t)k + 1 : 0;
    w.resize(nwin);
    if (nwin == 0) return;

    // 0=A, 1=C, 2=G, 3=T, 4=invalid
    static const uint8_t ENC[256] = {
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4, 4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,
        4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4, 4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    };

    // Precompute n * log2(n) for n = 0 .. k to avoid repeated log2 calls
    std::vector<double> nlgn((size_t)(k + 1));
    nlgn[0] = 0.0;
    for (int n = 1; n <= k; n++)
        nlgn[n] = (double)n * std::log2((double)n);

    const double inv_k      = 1.0 / (double)k;
    const double log2k      = std::log2((double)k);
    const double max_H      = 2.0;          // log2(4) = 2.0 for ACGT alphabet
    const double scale      = 1.0 / max_H;

    // Initialize counts for the first k-mer
    int cnt[5] = {0, 0, 0, 0, 0};
    for (int j = 0; j < k; j++)
        cnt[ENC[(uint8_t)seq[j]]]++;

    auto compute_weight = [&]() -> double {
        if (cnt[4] > 0) return w_min;   // contains N or invalid base
        double H = log2k - inv_k * (nlgn[cnt[0]] + nlgn[cnt[1]] +
                                     nlgn[cnt[2]] + nlgn[cnt[3]]);
        double weight = H * scale;
        return weight < w_min ? w_min : weight;
    };

    w[0] = compute_weight();
    for (uint64_t i = 1; i < nwin; i++) {
        cnt[ENC[(uint8_t)seq[i - 1]]]--;
        cnt[ENC[(uint8_t)seq[i + (uint64_t)k - 1]]]++;
        w[i] = compute_weight();
    }
}

#endif
