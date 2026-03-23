/**
 * eval_sketch_accuracy – synthetic benchmark for sketch Jaccard accuracy.
 *
 * Generates pairs of random DNA sequences with controlled per-base substitution
 * rates.  For each pair, builds all five sketch types and compares the estimated
 * Jaccard against the theoretical Mash formula:
 *
 *     J(p,k) = (1-p)^k / (2 - (1-p)^k)
 *
 * Sketch methods & their k-mer sizes:
 *   HyperLogLog  (k=32 hardcoded, 1024 registers)
 *   SetSketch    (k=32 hardcoded, 1024 registers)
 *   KSSD         (k=20, half_k=10 half_subk=6 drlevel=3)
 *   MinHash      (k=21, sketch size 1024)
 *   ProbMinHash4 (k=21, 1024 registers)
 *
 * §1  Fixed-length benchmark:  all pairs use seq_length (default 4 Mbp).
 * §2  Variable-length benchmark: lengths drawn log-uniformly from
 *     [500 Kbp, 20 Mbp], mimicking the broad size range of real genomes.
 *     Within each pair both sequences share the same length; length varies
 *     across pairs.
 *
 * Usage:
 *   exe_eval_sketch_acc [pairs_per_rate] [seq_length] [threads]
 *
 * Output:
 *   stdout  – two CSV sections (fixed-length then variable-length)
 *   stderr  – summary tables with MAE / RMSE / bias per method
 *
 * NOTE: Must be run from the examples/ directory so that the KSSD shuffle file
 *       (shuf_file/L3K10.shuf) is found.
 */

#include "Sketch.h"
#include "probmh.h"
#include "common.h"

#include <omp.h>
#include <sys/time.h>
#include <sys/stat.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

using namespace std;

static const char BASES[] = "ACGT";

static double theo_jaccard(double p, int k) {
    double q = pow(1.0 - p, k);
    return q / (2.0 - q);
}

static void gen_random_seq(char* buf, int len, mt19937_64& rng) {
    uniform_int_distribution<int> d4(0, 3);
    for (int i = 0; i < len; i++)
        buf[i] = BASES[d4(rng)];
    buf[len] = '\0';
}

static void gen_mutant(const char* base, char* out, int len,
                       double rate, mt19937_64& rng) {
    uniform_real_distribution<double> coin(0.0, 1.0);
    uniform_int_distribution<int>     d3(0, 2);
    for (int i = 0; i < len; i++) {
        if (coin(rng) < rate) {
            char alt[3]; int idx = 0;
            for (int b = 0; b < 4; b++)
                if (BASES[b] != base[i]) alt[idx++] = BASES[b];
            out[i] = alt[d3(rng)];
        } else {
            out[i] = base[i];
        }
    }
    out[len] = '\0';
}

struct PairResult {
    double rate;
    double theo_j32, theo_j21, theo_j20;
    double hll_j, ss_j, kssd_j, mh_j, pmh_j, oph_j, kmv_j;
    double tl1_j, tl2_j, tl4_j;   // Route C: Top-L truncated ProbMinHash
};

int main(int argc, char* argv[])
{
    int pairs_per_rate = (argc > 1) ? atoi(argv[1]) : 500;
    int seq_length     = (argc > 2) ? atoi(argv[2]) : 4000000;
    int numThreads     = (argc > 3) ? atoi(argv[3]) : 8;
    if (pairs_per_rate < 1) pairs_per_rate = 1;
    if (numThreads < 1) numThreads = 1;

    // ── KSSD shuffle file ──────────────────────────────────────────────────
    const char* shuf_path = "shuf_file/L3K10.shuf";
    struct stat st;
    if (stat(shuf_path, &st) != 0) {
        fprintf(stderr,
                "ERROR: KSSD shuffle file not found: %s\n"
                "       Run from the examples/ directory.\n", shuf_path);
        return 1;
    }
    Sketch::kssd_parameter_t kssdPara(10, 6, 3, shuf_path);

    vector<double> rates = {0.001, 0.005, 0.01, 0.02, 0.05,
                            0.1,   0.15,  0.2,  0.25, 0.3};
    int n_rates     = (int)rates.size();
    int total_pairs = n_rates * pairs_per_rate;

    fprintf(stderr,
            "=== Sketch Accuracy Evaluation ===\n"
            "  pairs/rate    : %d\n"
            "  seq length    : %d bp\n"
            "  mutation rates: %d levels (%.3f – %.3f)\n"
            "  total pairs   : %d  (%d sequences)\n"
            "  threads       : %d\n"
            "  sketch sizes  : HLL=1024  SS=1024  MH=1024  PMH=1024  OPH=1024  KMV=1024  TL{1,2,4}=1024\n"
            "  k-mer sizes   : HLL=32  SS=32  KSSD=20  MH=21  PMH=21  OPH=21  KMV=21  TL=21\n\n",
            pairs_per_rate, seq_length, n_rates,
            rates.front(), rates.back(),
            total_pairs, total_pairs * 2, numThreads);

    vector<PairResult> results(total_pairs);

    double t0 = get_sec();

    #pragma omp parallel for schedule(dynamic) num_threads(numThreads)
    for (int idx = 0; idx < total_pairs; idx++) {
        int  ri = idx / pairs_per_rate;
        double p = rates[ri];

        mt19937_64 rng(42ULL + (uint64_t)idx * 1000003ULL);

        vector<char> seq_a(seq_length + 1);
        vector<char> seq_b(seq_length + 1);
        gen_random_seq(seq_a.data(), seq_length, rng);
        gen_mutant(seq_a.data(), seq_b.data(), seq_length, p, rng);

        // ── HLL (k=32, np=10 → 1024 regs) ─────────────────────────────
        Sketch::HyperLogLog h1(10), h2(10);
        h1.update(seq_a.data());
        h2.update(seq_b.data());
        double hll_j = h1.jaccard_index(h2);

        // ── SetSketch (k=32, np=10 → 1024 regs) ───────────────────────
        Sketch::SetSketch s1(10), s2(10);
        s1.update(seq_a.data());
        s2.update(seq_b.data());
        double ss_j = s1.jaccard_index(s2);

        // ── KSSD (k=20, half_k=10, half_subk=6, drlevel=3) ────────────
        Sketch::Kssd* k1 = new Sketch::Kssd(kssdPara);
        Sketch::Kssd* k2 = new Sketch::Kssd(kssdPara);
        k1->update(seq_a.data());
        k2->update(seq_b.data());
        double kssd_j = k1->jaccard(k2);
        delete k1;
        delete k2;

        // ── MinHash (k=21, size=1024) ──────────────────────────────────
        Sketch::MinHash m1(21, 1024, 42, true), m2(21, 1024, 42, true);
        m1.update(seq_a.data());
        m2.update(seq_b.data());
        m1.finalize();
        m2.finalize();
        double mh_j = m1.jaccard(&m2);

        // ── ProbMinHash4 (m=1024, k=21) ───────────────────────────────
        Sketch::ProbMinHash4 pm1(1024, 21, 42), pm2(1024, 21, 42);
        pm1.update(seq_a.data(), seq_length);
        pm2.update(seq_b.data(), seq_length);
        double pmh_j = pm1.jaccard(pm2);

        // ── Route C: Top-L ProbMinHash (L=1, 2, 4) ─────────────────
        Sketch::ProbMinHash4 tl1a(1024,21,42,1), tl1b(1024,21,42,1);
        tl1a.update(seq_a.data(), seq_length);
        tl1b.update(seq_b.data(), seq_length);
        double tl1_j = tl1a.jaccard(tl1b);

        Sketch::ProbMinHash4 tl2a(1024,21,42,2), tl2b(1024,21,42,2);
        tl2a.update(seq_a.data(), seq_length);
        tl2b.update(seq_b.data(), seq_length);
        double tl2_j = tl2a.jaccard(tl2b);

        Sketch::ProbMinHash4 tl4a(1024,21,42,4), tl4b(1024,21,42,4);
        tl4a.update(seq_a.data(), seq_length);
        tl4b.update(seq_b.data(), seq_length);
        double tl4_j = tl4a.jaccard(tl4b);

        // ── ProbMinHash4OP – One-Permutation (m=1024, k=21) ─────────
        Sketch::ProbMinHash4OP op1(1024, 21, 42), op2(1024, 21, 42);
        op1.update(seq_a.data(), seq_length);
        op2.update(seq_b.data(), seq_length);
        double oph_j = op1.jaccard(op2);

        // ── ProbKMV – bottom-k (k=1024, kmer=21) ───────────────────
        Sketch::ProbKMV kv1(1024, 21, 42), kv2(1024, 21, 42);
        kv1.update(seq_a.data(), seq_length);
        kv2.update(seq_b.data(), seq_length);
        double kmv_j = kv1.jaccard(kv2);

        PairResult& r = results[idx];
        r.rate     = p;
        r.theo_j32 = theo_jaccard(p, 32);
        r.theo_j21 = theo_jaccard(p, 21);
        r.theo_j20 = theo_jaccard(p, 20);
        r.hll_j    = hll_j;
        r.ss_j     = ss_j;
        r.kssd_j   = kssd_j;
        r.mh_j     = mh_j;
        r.pmh_j    = pmh_j;
        r.oph_j    = oph_j;
        r.kmv_j    = kmv_j;
        r.tl1_j    = tl1_j;
        r.tl2_j    = tl2_j;
        r.tl4_j    = tl4_j;
    }

    double t1 = get_sec();
    fprintf(stderr, "Computation: %.2f s  (%.1f pairs/s)\n\n",
            t1 - t0, total_pairs / (t1 - t0));

    // ── CSV output ─────────────────────────────────────────────────────────
    printf("rate,theo_j_k32,theo_j_k21,theo_j_k20,"
           "hll_j,setsketch_j,kssd_j,minhash_j,probmh_j,oneperm_j,kmv_j,"
           "topl1_j,topl2_j,topl4_j\n");
    for (int i = 0; i < total_pairs; i++) {
        const PairResult& r = results[i];
        printf("%.4f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,"
               "%.6f,%.6f,%.6f\n",
               r.rate, r.theo_j32, r.theo_j21, r.theo_j20,
               r.hll_j, r.ss_j, r.kssd_j, r.mh_j, r.pmh_j, r.oph_j, r.kmv_j,
               r.tl1_j, r.tl2_j, r.tl4_j);
    }

    // ── Summary table ──────────────────────────────────────────────────────
    fprintf(stderr,
        "%-8s │ %-6s  %-33s │ %-6s  %-22s │ %-6s  %-33s %-22s\n",
        "rate",
        "thJ32", "HLL            SetSketch",
        "thJ20", "KSSD",
        "thJ21", "MinHash        ProbMH", "   OnePerm");
    fprintf(stderr,
        "%-8s │ %-6s  %-7s %-7s %-8s %-7s %-7s %-8s │ %-6s  %-7s %-7s %-8s │ %-6s  %-7s %-7s %-8s %-7s %-7s %-8s %-7s %-7s %-8s\n",
        "", "",
        "MAE", "RMSE", "bias", "MAE", "RMSE", "bias",
        "", "MAE", "RMSE", "bias",
        "", "MAE", "RMSE", "bias", "MAE", "RMSE", "bias",
        "MAE", "RMSE", "bias");
    fprintf(stderr, "─────────┼───────────────────────────────────────────"
                    "┼───────────────────────────────"
                    "┼──────────────────────────────────────────────"
                    "───────────────────────────\n");

    for (int ri = 0; ri < n_rates; ri++) {
        double p   = rates[ri];
        double j32 = theo_jaccard(p, 32);
        double j21 = theo_jaccard(p, 21);
        double j20 = theo_jaccard(p, 20);
        int    N   = pairs_per_rate;

        double hll_ae = 0, hll_se = 0, hll_bi = 0;
        double  ss_ae = 0,  ss_se = 0,  ss_bi = 0;
        double  kd_ae = 0,  kd_se = 0,  kd_bi = 0;
        double  mh_ae = 0,  mh_se = 0,  mh_bi = 0;
        double pmh_ae = 0, pmh_se = 0, pmh_bi = 0;
        double oph_ae = 0, oph_se = 0, oph_bi = 0;
        double kmv_ae = 0, kmv_se = 0, kmv_bi = 0;
        double tl1_ae = 0, tl1_se = 0, tl1_bi = 0;
        double tl2_ae = 0, tl2_se = 0, tl2_bi = 0;
        double tl4_ae = 0, tl4_se = 0, tl4_bi = 0;

        for (int pi = 0; pi < N; pi++) {
            const PairResult& r = results[ri * N + pi];

            double eh = r.hll_j - j32;
            hll_ae += fabs(eh); hll_se += eh * eh; hll_bi += eh;

            double es = r.ss_j - j32;
            ss_ae += fabs(es); ss_se += es * es; ss_bi += es;

            double ek = r.kssd_j - j20;
            kd_ae += fabs(ek); kd_se += ek * ek; kd_bi += ek;

            double em = r.mh_j - j21;
            mh_ae += fabs(em); mh_se += em * em; mh_bi += em;

            double ep = r.pmh_j - j21;
            pmh_ae += fabs(ep); pmh_se += ep * ep; pmh_bi += ep;

            double eo = r.oph_j - j21;
            oph_ae += fabs(eo); oph_se += eo * eo; oph_bi += eo;

            double ev = r.kmv_j - j21;
            kmv_ae += fabs(ev); kmv_se += ev * ev; kmv_bi += ev;

            double e1 = r.tl1_j - j21;
            tl1_ae += fabs(e1); tl1_se += e1 * e1; tl1_bi += e1;
            double e2 = r.tl2_j - j21;
            tl2_ae += fabs(e2); tl2_se += e2 * e2; tl2_bi += e2;
            double e4 = r.tl4_j - j21;
            tl4_ae += fabs(e4); tl4_se += e4 * e4; tl4_bi += e4;
        }

        fprintf(stderr,
            "%-8.3f │ %.4f  %.4f  %.4f  %+.4f  %.4f  %.4f  %+.4f"
            " │ %.4f  %.4f  %.4f  %+.4f"
            " │ %.4f  %.4f  %.4f  %+.4f  %.4f  %.4f  %+.4f  %.4f  %.4f  %+.4f  %.4f  %.4f  %+.4f\n",
            p, j32,
            hll_ae / N, sqrt(hll_se / N), hll_bi / N,
             ss_ae / N, sqrt( ss_se / N),  ss_bi / N,
            j20,
             kd_ae / N, sqrt( kd_se / N),  kd_bi / N,
            j21,
             mh_ae / N, sqrt( mh_se / N),  mh_bi / N,
            pmh_ae / N, sqrt(pmh_se / N), pmh_bi / N,
            oph_ae / N, sqrt(oph_se / N), oph_bi / N,
            kmv_ae / N, sqrt(kmv_se / N), kmv_bi / N);
    }

    // ── Global aggregates ──────────────────────────────────────────────────
    double g_hll = 0, g_ss = 0, g_kd = 0, g_mh = 0, g_pmh = 0, g_oph = 0, g_kmv = 0;
    double g_tl1 = 0, g_tl2 = 0, g_tl4 = 0;
    for (int i = 0; i < total_pairs; i++) {
        int ri = i / pairs_per_rate;
        double j32 = theo_jaccard(rates[ri], 32);
        double j21 = theo_jaccard(rates[ri], 21);
        double j20 = theo_jaccard(rates[ri], 20);
        g_hll += fabs(results[i].hll_j  - j32);
        g_ss  += fabs(results[i].ss_j   - j32);
        g_kd  += fabs(results[i].kssd_j - j20);
        g_mh  += fabs(results[i].mh_j   - j21);
        g_pmh += fabs(results[i].pmh_j  - j21);
        g_oph += fabs(results[i].oph_j  - j21);
        g_kmv += fabs(results[i].kmv_j  - j21);
        g_tl1 += fabs(results[i].tl1_j  - j21);
        g_tl2 += fabs(results[i].tl2_j  - j21);
        g_tl4 += fabs(results[i].tl4_j  - j21);
    }
    fprintf(stderr,
            "\nGlobal MAE:  HLL=%.5f  SetSketch=%.5f  KSSD=%.5f  MinHash=%.5f  ProbMH=%.5f  OnePerm=%.5f  KMV=%.5f\n",
            g_hll / total_pairs, g_ss / total_pairs, g_kd / total_pairs,
            g_mh / total_pairs, g_pmh / total_pairs, g_oph / total_pairs, g_kmv / total_pairs);

    // ── Route C: Top-L tradeoff table ───────────────────────────────────
    fprintf(stderr,
            "\n--- Route C: Top-L ProbMinHash Tradeoff (fixed-length) ---\n"
            "%-8s │ %8s  %8s  %8s  %8s\n",
            "rate", "ProbMH", "TL-4", "TL-2", "TL-1");
    fprintf(stderr,
            "─────────┼─────────────────────────────────────────\n");
    for (int ri = 0; ri < n_rates; ri++) {
        double p = rates[ri];
        int    N = pairs_per_rate;
        double a_pmh = 0, a_tl4 = 0, a_tl2 = 0, a_tl1 = 0;
        for (int pi = 0; pi < N; pi++) {
            const PairResult& r = results[ri * N + pi];
            double j = theo_jaccard(p, 21);
            a_pmh += fabs(r.pmh_j - j);
            a_tl4 += fabs(r.tl4_j - j);
            a_tl2 += fabs(r.tl2_j - j);
            a_tl1 += fabs(r.tl1_j - j);
        }
        fprintf(stderr, "%-8.3f │ %8.5f  %8.5f  %8.5f  %8.5f\n",
                p, a_pmh/N, a_tl4/N, a_tl2/N, a_tl1/N);
    }
    fprintf(stderr,
            "─────────┼─────────────────────────────────────────\n"
            "%-8s │ %8.5f  %8.5f  %8.5f  %8.5f\n",
            "GLOBAL", g_pmh/total_pairs, g_tl4/total_pairs,
            g_tl2/total_pairs, g_tl1/total_pairs);

    fprintf(stderr, "Fixed-len total time: %.2f s\n", get_sec() - t0);

    // ═══════════════════════════════════════════════════════════════════════
    // §2  VARIABLE-LENGTH BENCHMARK  (mixed-size pairs)
    //
    //     For each pair, seq_a and seq_b have *independently* sampled lengths
    //     L_a and L_b drawn log-uniformly from [500 Kbp, 20 Mbp].
    //     Let S = min(L_a, L_b).  seq_b is constructed as:
    //       • first S bases: substitution mutant of seq_a[0..S-1] at rate p
    //       • remaining |L_b − L_a| bases (if any): fresh random sequence
    //
    //     Ground-truth Jaccard:
    //       J(p, k, L_a, L_b) = S·(1−p)^k / (L_a + L_b − S·(1−p)^k)
    //     which degenerates to the standard formula when L_a = L_b.
    // ═══════════════════════════════════════════════════════════════════════

    struct VarPairResult {
        double rate;
        int    len_a, len_b;
        double theo_j32, theo_j21, theo_j20;
        double hll_j, ss_j, kssd_j, mh_j, pmh_j, oph_j, kmv_j;
        double tl1_j, tl2_j, tl4_j;
    };

    // Ground-truth Jaccard for unequal-length pairs
    auto theo_jaccard_vl = [](double p, int k, int La, int Lb) -> double {
        double S      = (double)min(La, Lb);
        double shared = S * pow(1.0 - p, k);
        return shared / (La + Lb - shared);
    };

    vector<VarPairResult> vr(total_pairs);

    // Pre-generate all (L_a, L_b) pairs – single-threaded for reproducibility
    {
        mt19937_64 lrng(0xDEADBEEFCAFEULL);
        uniform_real_distribution<double> log_ld(log(5e5), log(2e7));
        for (int i = 0; i < total_pairs; i++) {
            vr[i].len_a = (int)round(exp(log_ld(lrng)));
            vr[i].len_b = (int)round(exp(log_ld(lrng)));
        }
    }

    {
        int lo_a = vr[0].len_a, hi_a = vr[0].len_a;
        int lo_b = vr[0].len_b, hi_b = vr[0].len_b;
        double ratio_sum = 0;
        for (int i = 0; i < total_pairs; i++) {
            lo_a = min(lo_a, vr[i].len_a); hi_a = max(hi_a, vr[i].len_a);
            lo_b = min(lo_b, vr[i].len_b); hi_b = max(hi_b, vr[i].len_b);
            double ra = (double)vr[i].len_a, rb = (double)vr[i].len_b;
            ratio_sum += (ra > rb) ? ra / rb : rb / ra;
        }
        fprintf(stderr,
                "\n=== Variable-Length Benchmark (mixed-size pairs) ===\n"
                "  pairs/rate  : %d  (%d total)\n"
                "  L_a range   : %d – %d bp\n"
                "  L_b range   : %d – %d bp\n"
                "  mean |La/Lb|: %.2f\n"
                "  threads     : %d\n\n",
                pairs_per_rate, total_pairs,
                lo_a, hi_a, lo_b, hi_b,
                ratio_sum / total_pairs, numThreads);
    }

    double t2 = get_sec();

    #pragma omp parallel for schedule(dynamic) num_threads(numThreads)
    for (int idx = 0; idx < total_pairs; idx++) {
        int    ri2 = idx / pairs_per_rate;
        double p   = rates[ri2];
        int    La  = vr[idx].len_a;
        int    Lb  = vr[idx].len_b;
        int    S   = min(La, Lb);

        mt19937_64 rng(0xCAFEBABEULL + (uint64_t)idx * 1000033ULL);

        vector<char> seq_a(La + 1), seq_b(Lb + 1);

        // Build seq_a
        gen_random_seq(seq_a.data(), La, rng);

        // Build seq_b: first S bases = mutant of seq_a[0..S-1]
        gen_mutant(seq_a.data(), seq_b.data(), S, p, rng);
        // Tail extension (when Lb > La): fresh random bases
        if (Lb > S)
            gen_random_seq(seq_b.data() + S, Lb - S, rng);
        seq_b[Lb] = '\0';

        Sketch::HyperLogLog  h1(10), h2(10);
        h1.update(seq_a.data()); h2.update(seq_b.data());
        double hll_j = h1.jaccard_index(h2);

        Sketch::SetSketch    s1(10), s2(10);
        s1.update(seq_a.data()); s2.update(seq_b.data());
        double ss_j = s1.jaccard_index(s2);

        Sketch::Kssd* k1 = new Sketch::Kssd(kssdPara);
        Sketch::Kssd* k2 = new Sketch::Kssd(kssdPara);
        k1->update(seq_a.data()); k2->update(seq_b.data());
        double kssd_j = k1->jaccard(k2);
        delete k1; delete k2;

        Sketch::MinHash m1(21, 1024, 42, true), m2(21, 1024, 42, true);
        m1.update(seq_a.data()); m2.update(seq_b.data());
        m1.finalize(); m2.finalize();
        double mh_j = m1.jaccard(&m2);

        Sketch::ProbMinHash4 pm1(1024, 21, 42), pm2(1024, 21, 42);
        pm1.update(seq_a.data(), La); pm2.update(seq_b.data(), Lb);
        double pmh_j = pm1.jaccard(pm2);

        Sketch::ProbMinHash4 vtl1a(1024,21,42,1), vtl1b(1024,21,42,1);
        vtl1a.update(seq_a.data(), La); vtl1b.update(seq_b.data(), Lb);
        double vtl1_j = vtl1a.jaccard(vtl1b);

        Sketch::ProbMinHash4 vtl2a(1024,21,42,2), vtl2b(1024,21,42,2);
        vtl2a.update(seq_a.data(), La); vtl2b.update(seq_b.data(), Lb);
        double vtl2_j = vtl2a.jaccard(vtl2b);

        Sketch::ProbMinHash4 vtl4a(1024,21,42,4), vtl4b(1024,21,42,4);
        vtl4a.update(seq_a.data(), La); vtl4b.update(seq_b.data(), Lb);
        double vtl4_j = vtl4a.jaccard(vtl4b);

        Sketch::ProbMinHash4OP op1(1024, 21, 42), op2(1024, 21, 42);
        op1.update(seq_a.data(), La); op2.update(seq_b.data(), Lb);
        double oph_j = op1.jaccard(op2);

        Sketch::ProbKMV kv1(1024, 21, 42), kv2(1024, 21, 42);
        kv1.update(seq_a.data(), La); kv2.update(seq_b.data(), Lb);
        double kmv_j = kv1.jaccard(kv2);

        VarPairResult& r = vr[idx];
        r.rate     = p;
        r.theo_j32 = theo_jaccard_vl(p, 32, La, Lb);
        r.theo_j21 = theo_jaccard_vl(p, 21, La, Lb);
        r.theo_j20 = theo_jaccard_vl(p, 20, La, Lb);
        r.hll_j    = hll_j;
        r.ss_j     = ss_j;
        r.kssd_j   = kssd_j;
        r.mh_j     = mh_j;
        r.pmh_j    = pmh_j;
        r.oph_j    = oph_j;
        r.kmv_j    = kmv_j;
        r.tl1_j    = vtl1_j;
        r.tl2_j    = vtl2_j;
        r.tl4_j    = vtl4_j;
    }

    double t3 = get_sec();
    fprintf(stderr, "VarLen computation: %.2f s  (%.1f pairs/s)\n\n",
            t3 - t2, total_pairs / (t3 - t2));

    // ── Variable-length CSV output (stdout) ───────────────────────────────
    printf("\n# === VARIABLE-LENGTH RESULTS (mixed-size pairs) ===\n");
    printf("len_a_bp,len_b_bp,rate,theo_j_k32,theo_j_k21,theo_j_k20,"
           "hll_j,setsketch_j,kssd_j,minhash_j,probmh_j,oneperm_j,kmv_j,"
           "topl1_j,topl2_j,topl4_j\n");
    for (int i = 0; i < total_pairs; i++) {
        const VarPairResult& r = vr[i];
        printf("%d,%d,%.4f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,"
               "%.6f,%.6f,%.6f\n",
               r.len_a, r.len_b, r.rate,
               r.theo_j32, r.theo_j21, r.theo_j20,
               r.hll_j, r.ss_j, r.kssd_j, r.mh_j, r.pmh_j, r.oph_j, r.kmv_j,
               r.tl1_j, r.tl2_j, r.tl4_j);
    }

    // ── Global MAE summary – all lengths mixed (stderr) ───────────────────
    double v_hll = 0, v_ss = 0, v_kd = 0, v_mh = 0, v_pmh = 0, v_oph = 0, v_kmv = 0;
    double v_hll_se = 0, v_ss_se = 0, v_kd_se = 0, v_mh_se = 0, v_pmh_se = 0, v_oph_se = 0, v_kmv_se = 0;
    double v_tl1 = 0, v_tl2 = 0, v_tl4 = 0;
    double v_tl1_se = 0, v_tl2_se = 0, v_tl4_se = 0;
    for (int i = 0; i < total_pairs; i++) {
        const VarPairResult& r = vr[i];
        double eh = r.hll_j  - r.theo_j32; v_hll  += fabs(eh); v_hll_se  += eh*eh;
        double es = r.ss_j   - r.theo_j32; v_ss   += fabs(es); v_ss_se   += es*es;
        double ek = r.kssd_j - r.theo_j20; v_kd   += fabs(ek); v_kd_se   += ek*ek;
        double em = r.mh_j   - r.theo_j21; v_mh   += fabs(em); v_mh_se   += em*em;
        double ep = r.pmh_j  - r.theo_j21; v_pmh  += fabs(ep); v_pmh_se  += ep*ep;
        double eo = r.oph_j  - r.theo_j21; v_oph  += fabs(eo); v_oph_se  += eo*eo;
        double ev = r.kmv_j  - r.theo_j21; v_kmv  += fabs(ev); v_kmv_se  += ev*ev;
        double e1 = r.tl1_j  - r.theo_j21; v_tl1  += fabs(e1); v_tl1_se  += e1*e1;
        double e2 = r.tl2_j  - r.theo_j21; v_tl2  += fabs(e2); v_tl2_se  += e2*e2;
        double e4 = r.tl4_j  - r.theo_j21; v_tl4  += fabs(e4); v_tl4_se  += e4*e4;
    }
    double N = (double)total_pairs;
    fprintf(stderr,
        "\n--- Variable-Length Global MAE (all mutation rates, all sizes mixed) ---\n"
        "%-10s │ %8s %9s %8s %8s %8s %8s %8s\n",
        "metric", "HLL", "SetSketch", "KSSD", "MinHash", "ProbMH", "OnePerm", "KMV");
    fprintf(stderr,
        "───────────┼──────────────────────────────────────────────────────────────────────────\n");
    fprintf(stderr,
        "%-10s │ %8.5f %9.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n", "MAE",
        v_hll/N, v_ss/N, v_kd/N, v_mh/N, v_pmh/N, v_oph/N, v_kmv/N);
    fprintf(stderr,
        "%-10s │ %8.5f %9.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n", "RMSE",
        sqrt(v_hll_se/N), sqrt(v_ss_se/N), sqrt(v_kd_se/N),
        sqrt(v_mh_se/N), sqrt(v_pmh_se/N), sqrt(v_oph_se/N), sqrt(v_kmv_se/N));

    // ── Route C: Variable-length Top-L tradeoff table ───────────────────
    fprintf(stderr,
            "\n--- Route C: Top-L ProbMinHash Tradeoff (variable-length) ---\n"
            "%-10s │ %8s  %8s  %8s  %8s\n",
            "metric", "ProbMH", "TL-4", "TL-2", "TL-1");
    fprintf(stderr,
            "───────────┼─────────────────────────────────────────\n");
    fprintf(stderr,
            "%-10s │ %8.5f  %8.5f  %8.5f  %8.5f\n", "MAE",
            v_pmh/N, v_tl4/N, v_tl2/N, v_tl1/N);
    fprintf(stderr,
            "%-10s │ %8.5f  %8.5f  %8.5f  %8.5f\n", "RMSE",
            sqrt(v_pmh_se/N), sqrt(v_tl4_se/N), sqrt(v_tl2_se/N), sqrt(v_tl1_se/N));

    fprintf(stderr, "\nVarLen total time: %.2f s\n", get_sec() - t2);

    return 0;
}
