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
 * Usage:
 *   exe_eval_sketch_acc [pairs_per_rate] [seq_length] [threads]
 *
 * Output:
 *   stdout  – CSV with per-pair Jaccard estimates
 *   stderr  – summary table with MAE / RMSE / bias per mutation rate per method
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
    double hll_j, ss_j, kssd_j, mh_j, pmh_j;
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
            "  sketch sizes  : HLL=1024  SS=1024  MH=1024  PMH=1024\n"
            "  k-mer sizes   : HLL=32  SS=32  KSSD=20  MH=21  PMH=21\n\n",
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
        pm1.update(seq_a.data());
        pm2.update(seq_b.data());
        double pmh_j = pm1.jaccard(pm2);

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
    }

    double t1 = get_sec();
    fprintf(stderr, "Computation: %.2f s  (%.1f pairs/s)\n\n",
            t1 - t0, total_pairs / (t1 - t0));

    // ── CSV output ─────────────────────────────────────────────────────────
    printf("rate,theo_j_k32,theo_j_k21,theo_j_k20,"
           "hll_j,setsketch_j,kssd_j,minhash_j,probmh_j\n");
    for (int i = 0; i < total_pairs; i++) {
        const PairResult& r = results[i];
        printf("%.4f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
               r.rate, r.theo_j32, r.theo_j21, r.theo_j20,
               r.hll_j, r.ss_j, r.kssd_j, r.mh_j, r.pmh_j);
    }

    // ── Summary table ──────────────────────────────────────────────────────
    fprintf(stderr,
        "%-8s │ %-6s  %-33s │ %-6s  %-22s │ %-6s  %-33s\n",
        "rate",
        "thJ32", "HLL            SetSketch",
        "thJ20", "KSSD",
        "thJ21", "MinHash        ProbMH");
    fprintf(stderr,
        "%-8s │ %-6s  %-7s %-7s %-8s %-7s %-7s %-8s │ %-6s  %-7s %-7s %-8s │ %-6s  %-7s %-7s %-8s %-7s %-7s %-8s\n",
        "", "",
        "MAE", "RMSE", "bias", "MAE", "RMSE", "bias",
        "", "MAE", "RMSE", "bias",
        "", "MAE", "RMSE", "bias", "MAE", "RMSE", "bias");
    fprintf(stderr, "─────────┼───────────────────────────────────────────"
                    "┼───────────────────────────────"
                    "┼──────────────────────────────────────────────\n");

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
        }

        fprintf(stderr,
            "%-8.3f │ %.4f  %.4f  %.4f  %+.4f  %.4f  %.4f  %+.4f"
            " │ %.4f  %.4f  %.4f  %+.4f"
            " │ %.4f  %.4f  %.4f  %+.4f  %.4f  %.4f  %+.4f\n",
            p, j32,
            hll_ae / N, sqrt(hll_se / N), hll_bi / N,
             ss_ae / N, sqrt( ss_se / N),  ss_bi / N,
            j20,
             kd_ae / N, sqrt( kd_se / N),  kd_bi / N,
            j21,
             mh_ae / N, sqrt( mh_se / N),  mh_bi / N,
            pmh_ae / N, sqrt(pmh_se / N), pmh_bi / N);
    }

    // ── Global aggregates ──────────────────────────────────────────────────
    double g_hll = 0, g_ss = 0, g_kd = 0, g_mh = 0, g_pmh = 0;
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
    }
    fprintf(stderr,
            "\nGlobal MAE:  HLL=%.5f  SetSketch=%.5f  KSSD=%.5f  MinHash=%.5f  ProbMH=%.5f\n",
            g_hll / total_pairs, g_ss / total_pairs, g_kd / total_pairs,
            g_mh / total_pairs, g_pmh / total_pairs);
    fprintf(stderr, "Total time: %.2f s\n", get_sec() - t0);

    return 0;
}
