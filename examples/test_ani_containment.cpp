/**
 * test_ani_containment – verify ANI and containment for all four algorithms.
 *
 * Three test scenarios:
 *   (A) Equal-size sequences A and B, related by mutation rate p.
 *       Expected: ANI ≈ 1-p,  C(A⊆B) ≈ C(B⊆A) ≈ (1-p)^k = 2J/(1+J)
 *   (B) A is a prefix of C = A ‖ random_suffix  (A entirely contained in C).
 *       Expected: C(A⊆C) ≈ 1.0,  C(C⊆A) ≈ 0.5
 *   (C) Identical sequences.
 *       Expected: ANI = 1.0,  containment = 1.0,  Jaccard = 1.0
 */

#include "Sketch.h"
#include "fastkmv.h"
#include "probmh.h"
#include "common.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <vector>

using namespace std;

static const char BASES[] = "ACGT";

static void gen_seq(char* buf, int len, mt19937_64& rng) {
    uniform_int_distribution<int> d4(0, 3);
    for (int i = 0; i < len; i++) buf[i] = BASES[d4(rng)];
    buf[len] = '\0';
}

static void gen_mutant(const char* src, char* dst, int len,
                       double rate, mt19937_64& rng) {
    uniform_real_distribution<double> coin(0.0, 1.0);
    uniform_int_distribution<int>     d3(0, 2);
    for (int i = 0; i < len; i++) {
        if (coin(rng) < rate) {
            int orig = (src[i]=='A'?0 : src[i]=='C'?1 : src[i]=='G'?2 : 3);
            dst[i] = BASES[(orig + 1 + d3(rng)) % 4];
        } else {
            dst[i] = src[i];
        }
    }
    dst[len] = '\0';
}

static double theo_J(double p, int k) {
    double q = pow(1.0 - p, k);
    return q / (2.0 - q);
}

// for equal-size sequences: C(A⊆B) = 2J/(1+J) = (1-p)^k
static double theo_C_equal(double p, int k) {
    return pow(1.0 - p, k);
}

// ── print one row ────────────────────────────────────────────────────────────
static void row(const char* label, double got, double expected,
                double tol = 0.05) {
    bool ok = fabs(got - expected) <= tol;
    printf("  %-30s got=%.5f  expected=%.5f  %s\n",
           label, got, expected, ok ? "[OK]" : "[WARN]");
}

// ── main ─────────────────────────────────────────────────────────────────────
int main() {
    const int    L   = 1000000;   // 1 M bp
    const double P   = 0.05;      // 5 % substitution rate
    const int    K21 = 21;
    const int    K32 = 32;

    mt19937_64 rng(42);

    vector<char> seqA(L + 1), seqB(L + 1), seqC(2*L + 1), seqD(L + 1);
    gen_seq(seqA.data(), L, rng);
    gen_mutant(seqA.data(), seqB.data(), L, P, rng);
    // seqC = seqA ‖ random suffix  (A entirely contained in C)
    memcpy(seqC.data(), seqA.data(), L);
    gen_seq(seqC.data() + L, L, rng);
    seqC[2*L] = '\0';
    // seqD = seqA (identical copy)
    memcpy(seqD.data(), seqA.data(), L + 1);

    const double tj21 = theo_J(P, K21);
    const double tj32 = theo_J(P, K32);
    const double tc21 = theo_C_equal(P, K21);
    const double tc32 = theo_C_equal(P, K32);
    const double tani = 1.0 - P;

    printf("=== Theoretical reference values (mut_rate=%.2f, L=%d bp) ===\n",
           P, L);
    printf("  true ANI                = %.6f\n", tani);
    printf("  theo Jaccard  k=21      = %.6f\n", tj21);
    printf("  theo Jaccard  k=32      = %.6f\n", tj32);
    printf("  theo C(A⊆B) equal k=21 = %.6f  (= (1-p)^k = 2J/(1+J))\n", tc21);
    printf("  theo C(A⊆B) equal k=32 = %.6f\n", tc32);
    printf("\n");

    // ════════════════════════════════════════════════════════════════════════
    printf("────────────── FastKMV (k=1024, kmer=21) ──────────────\n");
    {
        Sketch::FastKMV kvA(1024, K21, 42), kvB(1024, K21, 42);
        Sketch::FastKMV kvC(1024, K21, 42), kvD(1024, K21, 42);
        kvA.update(seqA.data(), L);
        kvB.update(seqB.data(), L);
        kvC.update(seqC.data(), 2*L);
        kvD.update(seqD.data(), L);

        printf("  [A] Equal-size mutated pair\n");
        row("Jaccard(A,B)",        kvA.jaccard(kvB),        tj21);
        row("ANI(A,B)",            kvA.ani(kvB),            tani);
        row("C(A⊆B) equal-size",  kvA.containment(kvB),    tc21);
        row("C(B⊆A) equal-size",  kvB.containment(kvA),    tc21);

        printf("  [B] Subset: A is prefix of C=A‖random\n");
        row("C(A⊆C) [expect ~1.0]",kvA.containment(kvC),   1.0, 0.10);
        row("C(C⊆A) [expect ~0.5]",kvC.containment(kvA),   0.5, 0.10);

        printf("  [C] Identical sequences\n");
        row("Jaccard(A,D)",        kvA.jaccard(kvD),        1.0, 0.01);
        row("ANI(A,D)",            kvA.ani(kvD),            1.0, 0.01);
        row("C(A⊆D) identical",   kvA.containment(kvD),    1.0, 0.01);

        printf("  Cardinality check\n");
        row("card(A)  [true~1M]",  kvA.cardinality(), (double)(L-K21+1), 0.05*(L-K21+1));
        row("card(C)  [true~2M]",  kvC.cardinality(), (double)(2*L-K21+1), 0.05*(2*L-K21+1));
    }
    printf("\n");

    // ════════════════════════════════════════════════════════════════════════
    printf("────────────── HLL (np=13, kmer=32) ──────────────\n");
    {
        Sketch::HyperLogLog hA(13), hB(13), hC(13), hD(13);
        hA.update(seqA.data());
        hB.update(seqB.data());
        hC.update(seqC.data());
        hD.update(seqD.data());

        printf("  [A] Equal-size mutated pair\n");
        row("Jaccard(A,B)",        hA.jaccard_index(hB),    tj32);
        row("ANI(A,B)",            hA.ani(hB, K32),         tani);
        row("C(A⊆B) equal-size",  hA.containment(hB),      tc32);
        row("C(B⊆A) equal-size",  hB.containment(hA),      tc32);

        printf("  [B] Subset: A is prefix of C=A‖random\n");
        row("C(A⊆C) [expect ~1.0]",hA.containment(hC),     1.0, 0.10);
        row("C(C⊆A) [expect ~0.5]",hC.containment(hA),     0.5, 0.10);

        printf("  [C] Identical sequences\n");
        row("Jaccard(A,D)",        hA.jaccard_index(hD),    1.0, 0.01);
        row("ANI(A,D)",            hA.ani(hD, K32),         1.0, 0.01);
        row("C(A⊆D) identical",   hA.containment(hD),      1.0, 0.01);

        printf("  Cardinality check\n");
        row("card(A)  [true~1M]",  hA.cardinality(), (double)(L-K32+1), 0.05*(L-K32+1));
        row("card(C)  [true~2M]",  hC.cardinality(), (double)(2*L-K32+1), 0.05*(2*L-K32+1));
    }
    printf("\n");

    // ════════════════════════════════════════════════════════════════════════
    printf("────────────── SetSketch (np=13, kmer=32) ──────────────\n");
    {
        Sketch::SetSketch sA(13), sB(13), sC(13), sD(13);
        sA.update(seqA.data());
        sB.update(seqB.data());
        sC.update(seqC.data());
        sD.update(seqD.data());

        printf("  [A] Equal-size mutated pair\n");
        row("Jaccard(A,B)",        sA.jaccard_index(sB),    tj32);
        row("ANI(A,B)",            sA.ani(sB, K32),         tani);
        row("C(A⊆B) equal-size",  sA.containment(sB),      tc32);
        row("C(B⊆A) equal-size",  sB.containment(sA),      tc32);

        printf("  [B] Subset: A is prefix of C=A‖random\n");
        row("C(A⊆C) [expect ~1.0]",sA.containment(sC),     1.0, 0.10);
        row("C(C⊆A) [expect ~0.5]",sC.containment(sA),     0.5, 0.10);

        printf("  [C] Identical sequences\n");
        row("Jaccard(A,D)",        sA.jaccard_index(sD),    1.0, 0.01);
        row("ANI(A,D)",            sA.ani(sD, K32),         1.0, 0.01);
        row("C(A⊆D) identical",   sA.containment(sD),      1.0, 0.01);

        printf("  Cardinality check\n");
        row("card(A)  [true~1M]",  sA.cardinality(), (double)(L-K32+1), 0.05*(L-K32+1));
        row("card(C)  [true~2M]",  sC.cardinality(), (double)(2*L-K32+1), 0.05*(2*L-K32+1));
    }
    printf("\n");

    // ════════════════════════════════════════════════════════════════════════
    printf("────────────── ProbMinHash4 (m=1024, kmer=21) ──────────────\n");
    {
        Sketch::ProbMinHash4 pA(1024, K21, 42), pB(1024, K21, 42);
        Sketch::ProbMinHash4 pC(1024, K21, 42), pD(1024, K21, 42);
        pA.update(seqA.data(), L);
        pB.update(seqB.data(), L);
        pC.update(seqC.data(), 2*L);
        pD.update(seqD.data(), L);

        printf("  [A] Equal-size mutated pair\n");
        row("Jaccard(A,B)",        pA.jaccard(pB),          tj21);
        row("ANI(A,B)",            pA.ani(pB),              tani);
        row("C(A⊆B) equal-size",  pA.containment(pB),      tc21);
        row("C(B⊆A) equal-size",  pB.containment(pA),      tc21);

        printf("  [B] Subset: A is prefix of C=A‖random\n");
        row("C(A⊆C) [expect ~1.0]",pA.containment(pC),     1.0, 0.10);
        row("C(C⊆A) [expect ~0.5]",pC.containment(pA),     0.5, 0.10);

        printf("  [C] Identical sequences\n");
        row("Jaccard(A,D)",        pA.jaccard(pD),          1.0, 0.01);
        row("ANI(A,D)",            pA.ani(pD),              1.0, 0.01);
        row("C(A⊆D) identical",   pA.containment(pD),      1.0, 0.01);

        printf("  total_weight check\n");
        row("total_weight(A) [~1M]", pA.total_weight(),
            (double)(L - K21 + 1), 0.01*(L - K21 + 1));
        row("total_weight(C) [~2M]", pC.total_weight(),
            (double)(2*L - K21 + 1), 0.01*(2*L - K21 + 1));
    }
    printf("\n");

    return 0;
}
