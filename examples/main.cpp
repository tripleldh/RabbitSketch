/**
 * rabbitsketch – unified CLI front-end for the RabbitSketch library.
 *
 * Algorithms (mutually exclusive):
 *   --minhash      classical bottom-k MinHash
 *   --kssd         dimensionality-reduced k-mer sketch
 *   --hll          HyperLogLog cardinality / Jaccard sketch
 *   --probminhash  weighted ProbMinHash4 (entropy-weighted update)
 *   --bindash      b-bit one-permutation hashing
 *   --setsketch    register-based set sketch
 *   --fastkmv      FastKMV bottom-k sketch
 *
 * Modes:
 *   pair mode      -i fileA fileB                 (two FASTA/FASTQ inputs)
 *   list mode      -i list.txt -l                 (file list, all-pairs)
 *
 *   --index        in list mode, run the inverted-index pipeline
 *                  (currently supported for: --fastkmv --probminhash
 *                                            --kssd     --setsketch)
 *
 * Common options:
 *   -o <path>      output file (default: rabbitsketch.dist)
 *   -k <int>       k-mer size (default per algorithm)
 *   -t <int>       OpenMP threads     (default: omp_get_max_threads())
 *   -d <float>     max Mash distance for filter (default 0.05; use 1.0 for no filter)
 *   --seed <u64>   hash seed          (default 42)
 *
 * Algorithm-specific options:
 *   --minhash       -s   <int>     sketch size                    (default 1000)
 *   --kssd          --drlevel <int>                                (default 3)
 *   --hll           -b   <int>     log2(register count)            (default 13)
 *   --probminhash   -m   <int>     registers                       (default 1024)
 *                   -L   <int>     max-L truncation (0 = full)     (default 0)
 *                   --no-entropy   uniform-weighted update
 *                                  (entropy mode is on by default,
 *                                   matching test_ProbMinHash.cpp)
 *   --bindash       --bbits <int>  bits per bin                    (default 16)
 *                   --sketch64 <int> 64-bin groups                 (default 32)
 *   --setsketch     -b   <int>     log2(register count)            (default 13)
 *                   -a   <float>   parameter a                     (default 20.0)
 *                   -B   <float>   base                            (default 2.0)
 *   --fastkmv       -K   <int>     sketch size (k smallest)        (default 1024)
 *
 * Output (TSV):  <fileA>\t<fileB>\t<distance>
 * Distance is the Mash distance D = -ln(2J/(1+J))/k for all algorithms.
 */

#include "Sketch.h"
#include "BinDash.h"
#include "fastkmv.h"
#include "probmh.h"
#include "SetSketch.h"
#include "InvertedIndex.h"
#include "common.h"
#include "kseq.h"
#include "phmap.h"

#include <zlib.h>
#include <omp.h>
#include <err.h>
#include <sys/stat.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <tuple>
#include <numeric>
#include <string>
#include <vector>

KSEQ_INIT(gzFile, gzread)

// ─────────────────────────────────────────────────────────────────────────────
//  Argument parsing
// ─────────────────────────────────────────────────────────────────────────────
enum class Algo { NONE, MINHASH, KSSD, HLL, PROBMINHASH, BINDASH, SETSKETCH, FASTKMV };

struct Args {
    Algo        algo        = Algo::NONE;
    bool        useIndex    = false;
    bool        listMode    = false;
    std::vector<std::string> inputs;
    std::string output      = "rabbitsketch.dist";

    int         kmerSize    = 21;
    int         threads     = omp_get_max_threads();
    double      maxDist     = 0.05;  // typical genomic Mash-distance cutoff
    uint64_t    seed        = 42;

    int         minhashSize = 1000;
    int         kssdDrlevel = 3;
    int         hllBits     = 13;

    uint32_t    pmhM        = 1024;
    uint32_t    pmhMaxL     = 0;
    bool        pmhEntropy  = true;   // matches test_ProbMinHash.cpp default

    uint32_t    bdBbits     = 16;
    uint32_t    bdSketch64  = 32;

    int         ssBits      = 13;
    double      ssA         = 20.0;
    double      ssBase      = 2.0;

    uint32_t    fkmvK       = 1024;
};

static void print_usage(const char* prog) {
    std::cerr <<
"Usage: " << prog << " --<algo> [--index] -i <file>... [-l] [options]\n\n"
"Algorithms (choose one):\n"
"  --minhash --kssd --hll --probminhash --bindash --setsketch --fastkmv\n\n"
"Modes:\n"
"  -i fileA fileB         pairwise distance between two FASTA/FASTQ files\n"
"  -i list.txt -l         all-pairs over the file list\n"
"  --index                use inverted-index pipeline (list mode only;\n"
"                         supported for fastkmv/probminhash/kssd/setsketch)\n\n"
"Common options:\n"
"  -o <path>              output file (default: rabbitsketch.dist)\n"
"  -k <int>               k-mer size (default 21; setsketch uses 32)\n"
"  -t <int>               threads (default: omp_get_max_threads())\n"
"  -d <float>             max Mash distance filter (default 0.05; use 1.0 for no filter)\n"
"  --seed <u64>           hash seed (default 42)\n\n"
"Algorithm-specific:\n"
"  --minhash      -s <int>             sketch size      (1000)\n"
"  --kssd         --drlevel <int>      (3)\n"
"  --hll          -b <int>             log2(registers)  (13)\n"
"  --probminhash  -m <int> -L <int> [--no-entropy] (1024, 0; entropy on by default)\n"
"  --bindash      --bbits <int> --sketch64 <int>  (16, 32)\n"
"  --setsketch    -b <int> -a <float> -B <float>  (13, 20.0, 2.0)\n"
"  --fastkmv      -K <int>             sketch size      (1024)\n";
}

static int parse_args(int argc, char* argv[], Args& a) {
    auto need = [&](int i, const char* opt) {
        if (i + 1 >= argc) {
            std::cerr << "ERROR: option " << opt << " requires an argument\n";
            std::exit(1);
        }
    };
    auto setAlgo = [&](Algo x) {
        if (a.algo != Algo::NONE) {
            std::cerr << "ERROR: more than one --<algo> flag specified\n";
            std::exit(1);
        }
        a.algo = x;
    };

    for (int i = 1; i < argc; ++i) {
        std::string s = argv[i];
        if      (s == "--minhash")      setAlgo(Algo::MINHASH);
        else if (s == "--kssd")         setAlgo(Algo::KSSD);
        else if (s == "--hll")          setAlgo(Algo::HLL);
        else if (s == "--probminhash")  setAlgo(Algo::PROBMINHASH);
        else if (s == "--bindash")      setAlgo(Algo::BINDASH);
        else if (s == "--setsketch")  { setAlgo(Algo::SETSKETCH); a.kmerSize = 32; }
        else if (s == "--fastkmv")      setAlgo(Algo::FASTKMV);
        else if (s == "--index")        a.useIndex = true;
        else if (s == "-l")             a.listMode = true;
        else if (s == "-i") {
            need(i, "-i");
            const size_t before = a.inputs.size();
            while (i + 1 < argc && argv[i + 1][0] != '-')
                a.inputs.emplace_back(argv[++i]);
            if (a.inputs.size() == before) {
                std::cerr << "ERROR: -i must be followed by at least one file "
                             "(got next token: '" << argv[i + 1] << "')\n";
                std::exit(1);
            }
        }
        else if (s == "-o")             { need(i, "-o"); a.output    = argv[++i]; }
        else if (s == "-k")             { need(i, "-k"); a.kmerSize  = std::stoi(argv[++i]); }
        else if (s == "-t")             { need(i, "-t"); a.threads   = std::max(1, std::stoi(argv[++i])); }
        else if (s == "-d")             { need(i, "-d"); a.maxDist   = std::stod(argv[++i]); }
        else if (s == "--seed")         { need(i, "--seed"); a.seed  = std::stoull(argv[++i]); }
        else if (s == "-s")             { need(i, "-s"); a.minhashSize = std::stoi(argv[++i]); }
        else if (s == "--drlevel")      { need(i, "--drlevel"); a.kssdDrlevel = std::stoi(argv[++i]); }
        else if (s == "-b")             { need(i, "-b"); a.hllBits = a.ssBits = std::stoi(argv[++i]); }
        else if (s == "-m")             { need(i, "-m"); a.pmhM = static_cast<uint32_t>(std::stoul(argv[++i])); }
        else if (s == "-L")             { need(i, "-L"); a.pmhMaxL = static_cast<uint32_t>(std::stoul(argv[++i])); }
        else if (s == "--entropy")      a.pmhEntropy = true;
        else if (s == "--no-entropy")   a.pmhEntropy = false;
        else if (s == "--bbits")        { need(i, "--bbits");    a.bdBbits    = static_cast<uint32_t>(std::stoul(argv[++i])); }
        else if (s == "--sketch64")     { need(i, "--sketch64"); a.bdSketch64 = static_cast<uint32_t>(std::stoul(argv[++i])); }
        else if (s == "-a")             { need(i, "-a"); a.ssA    = std::stod(argv[++i]); }
        else if (s == "-B")             { need(i, "-B"); a.ssBase = std::stod(argv[++i]); }
        else if (s == "-K")             { need(i, "-K"); a.fkmvK  = static_cast<uint32_t>(std::stoul(argv[++i])); }
        else if (s == "-h" || s == "--help") { print_usage(argv[0]); std::exit(0); }
        else {
            std::cerr << "ERROR: unknown option: " << s << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    if (a.algo == Algo::NONE) {
        std::cerr << "ERROR: no --<algo> selected\n";
        print_usage(argv[0]);
        return 1;
    }
    if (a.inputs.empty()) {
        std::cerr << "ERROR: -i is required\n";
        return 1;
    }
    if (a.listMode) {
        if (a.inputs.size() != 1) {
            std::cerr << "ERROR: -l mode expects exactly one file (the list)\n";
            return 1;
        }
    } else {
        if (a.inputs.size() != 2) {
            std::cerr << "ERROR: pair mode requires exactly two files (or use -l)\n";
            return 1;
        }
    }
    if (a.useIndex && !a.listMode) {
        std::cerr << "WARNING: --index is only meaningful in list (-l) mode; ignored\n";
        a.useIndex = false;
    }
    return 0;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Helpers
// ─────────────────────────────────────────────────────────────────────────────
static std::vector<std::string> load_file_list(const std::string& path) {
    std::ifstream fs(path);
    if (!fs) err(errno, "cannot open file list: %s", path.c_str());
    std::vector<std::string> out;
    std::string line;
    while (std::getline(fs, line))
        if (!line.empty()) out.push_back(line);
    return out;
}

static inline double mash_distance(double jaccard, int kmer) {
    if (jaccard >= 1.0) return 0.0;
    if (jaccard <= 0.0) return 1.0;
    return -std::log(2.0 * jaccard / (1.0 + jaccard)) / static_cast<double>(kmer);
}

static void write_pair(const std::string& outPath,
                       const std::string& a, const std::string& b, double d)
{
    FILE* fp = std::fopen(outPath.c_str(), "w");
    if (!fp) err(errno, "cannot open output: %s", outPath.c_str());
    std::fprintf(fp, "%s\t%s\t%.6f\n", a.c_str(), b.c_str(), d);
    std::fclose(fp);
    std::cerr << "output: " << outPath << "\n";
}

// Stream a FASTA/FASTQ file and invoke `cb(seq, len)` for every record.
template<class CB>
static void stream_seq(const std::string& path, CB&& cb) {
    gzFile fp = gzopen(path.c_str(), "r");
    if (!fp) { std::cerr << "WARN: cannot open " << path << "\n"; return; }
    kseq_t* ks = kseq_init(fp);
    while (kseq_read(ks) >= 0)
        cb(ks->seq.s, static_cast<uint64_t>(ks->seq.l));
    kseq_destroy(ks);
    gzclose(fp);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Per-algorithm sketch builders – all return heap-allocated T*
// ─────────────────────────────────────────────────────────────────────────────
static Sketch::MinHash* build_minhash(const std::string& path, const Args& a) {
    auto* sk = new Sketch::MinHash(a.kmerSize, a.minhashSize,
                                   static_cast<uint32_t>(a.seed), true);
    sk->fileName = path;
    stream_seq(path, [&](char* seq, uint64_t /*len*/) { sk->update(seq); });
    sk->finalize();
    return sk;
}

static Sketch::HyperLogLog* build_hll(const std::string& path, const Args& a) {
    auto* sk = new Sketch::HyperLogLog(a.hllBits);
    stream_seq(path, [&](char* seq, uint64_t /*len*/) { sk->update(seq); });
    return sk;
}

// Build (or cache) a Kssd parameter object.  generate_shuffle_dim() uses the
// global libc PRNG (srand/rand), which is *not* thread-safe; if many builders
// construct kssd_parameter_t in parallel the shuffles collide and the dim map
// ends up partially empty, producing zero-Jaccard sketches.  We therefore
// build one parameter object per (half_k, half_subk, drlevel) tuple under a
// mutex and reuse it across all threads.
static const Sketch::kssd_parameter_t&
get_kssd_params(int half_k, int half_subk, int drlevel) {
    static std::mutex mu;
    static std::vector<std::pair<std::tuple<int,int,int>,
                                 std::shared_ptr<Sketch::kssd_parameter_t>>> cache;
    std::lock_guard<std::mutex> lk(mu);
    const auto key = std::make_tuple(half_k, half_subk, drlevel);
    for (auto& [k, v] : cache) if (k == key) return *v;
    cache.emplace_back(key,
        std::make_shared<Sketch::kssd_parameter_t>(half_k, half_subk, drlevel));
    return *cache.back().second;
}

static Sketch::Kssd* build_kssd(const std::string& path, const Args& a) {
    auto* sk = new Sketch::Kssd(
        get_kssd_params(/*half_k*/ a.kmerSize / 2,
                        /*half_subk*/ 6,
                        /*drlevel*/ a.kssdDrlevel));
    sk->fileName = path;
    // Kssd::update internally calls SetToList() at the end of every call, so
    // after the loop hashList is already merged + sorted and ready for
    // jaccard()/distance().
    stream_seq(path, [&](char* seq, uint64_t /*len*/) { sk->update(seq); });
    return sk;
}

static Sketch::ProbMinHash4* build_probminhash(const std::string& path, const Args& a) {
    auto* sk = new Sketch::ProbMinHash4(a.pmhM, a.kmerSize, a.seed, a.pmhMaxL);
    stream_seq(path, [&](char* seq, uint64_t len) {
        if (len < static_cast<uint64_t>(a.kmerSize)) return;
        if (a.pmhEntropy) sk->updateEntropy(seq, len);
        else              sk->update(seq, len);
    });
    return sk;
}

static Sketch::BinDash* build_bindash(const std::string& path, const Args& a) {
    auto* sk = new Sketch::BinDash(a.bdSketch64, a.kmerSize, a.bdBbits, a.seed);
    stream_seq(path, [&](char* seq, uint64_t len) { sk->update(seq, len); });
    sk->finalize();
    return sk;
}

static Sketch::SetSketch* build_setsketch(const std::string& path, const Args& a) {
    auto* sk = new Sketch::SetSketch(a.ssBits, a.ssBase, a.ssA);
    stream_seq(path, [&](char* seq, uint64_t len) {
        sk->update(seq, static_cast<size_t>(len));
    });
    return sk;
}

static Sketch::FastKMV* build_fastkmv(const std::string& path, const Args& a) {
    auto* sk = new Sketch::FastKMV(a.fkmvK, a.kmerSize, a.seed);
    stream_seq(path, [&](char* seq, uint64_t len) { sk->update(seq, len); });
    return sk;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Pair-mode: build two sketches, write one distance line
// ─────────────────────────────────────────────────────────────────────────────
static int run_pair(const Args& a) {
    const std::string& fA = a.inputs[0];
    const std::string& fB = a.inputs[1];
    double d = 1.0;

    switch (a.algo) {
    case Algo::MINHASH: {
        std::unique_ptr<Sketch::MinHash> sa(build_minhash(fA, a));
        std::unique_ptr<Sketch::MinHash> sb(build_minhash(fB, a));
        d = sa->distance(sb.get());
        break;
    }
    case Algo::HLL: {
        std::unique_ptr<Sketch::HyperLogLog> sa(build_hll(fA, a));
        std::unique_ptr<Sketch::HyperLogLog> sb(build_hll(fB, a));
        d = mash_distance(sa->jaccard_index(*sb), a.kmerSize);
        break;
    }
    case Algo::KSSD: {
        std::unique_ptr<Sketch::Kssd> sa(build_kssd(fA, a));
        std::unique_ptr<Sketch::Kssd> sb(build_kssd(fB, a));
        d = sa->distance(sb.get());
        break;
    }
    case Algo::PROBMINHASH: {
        std::unique_ptr<Sketch::ProbMinHash4> sa(build_probminhash(fA, a));
        std::unique_ptr<Sketch::ProbMinHash4> sb(build_probminhash(fB, a));
        d = mash_distance(sa->jaccard(*sb), a.kmerSize);
        break;
    }
    case Algo::BINDASH: {
        std::unique_ptr<Sketch::BinDash> sa(build_bindash(fA, a));
        std::unique_ptr<Sketch::BinDash> sb(build_bindash(fB, a));
        d = mash_distance(sa->jaccard(*sb), a.kmerSize);
        break;
    }
    case Algo::SETSKETCH: {
        std::unique_ptr<Sketch::SetSketch> sa(build_setsketch(fA, a));
        std::unique_ptr<Sketch::SetSketch> sb(build_setsketch(fB, a));
        d = mash_distance(sa->jaccard_index(*sb), a.kmerSize);
        break;
    }
    case Algo::FASTKMV: {
        std::unique_ptr<Sketch::FastKMV> sa(build_fastkmv(fA, a));
        std::unique_ptr<Sketch::FastKMV> sb(build_fastkmv(fB, a));
        d = mash_distance(sa->jaccard(*sb), a.kmerSize);
        break;
    }
    default: return 1;
    }

    write_pair(a.output, fA, fB, d);
    return 0;
}

// ─────────────────────────────────────────────────────────────────────────────
//  List-mode: all-pairs (no inverted index)
// ─────────────────────────────────────────────────────────────────────────────
template<class T, class BuildFn, class DistFn>
static void list_allpairs(const Args& a,
                          const std::vector<std::string>& files,
                          BuildFn build, DistFn dist)
{
    const int N = static_cast<int>(files.size());
    std::vector<std::unique_ptr<T>> sks(N);

    double t0 = get_sec();
    #pragma omp parallel for num_threads(a.threads) schedule(dynamic)
    for (int i = 0; i < N; ++i) sks[i].reset(build(files[i], a));
    std::cerr << "sketch time: " << get_sec() - t0 << " s\n";

    double t1 = get_sec();
    std::vector<std::string> bufs(a.threads);

    #pragma omp parallel for num_threads(a.threads) schedule(dynamic, 1)
    for (int i = 0; i < N; ++i) {
        const int tid = omp_get_thread_num();
        std::string& buf = bufs[tid];
        for (int j = i + 1; j < N; ++j) {
            const double d = dist(sks[i].get(), sks[j].get());
            if (d < a.maxDist) {
                char line[2048];
                int len = std::snprintf(line, sizeof(line), "%s\t%s\t%.6f\n",
                                        files[i].c_str(), files[j].c_str(), d);
                buf.append(line, static_cast<size_t>(len));
            }
        }
    }

    FILE* fp = std::fopen(a.output.c_str(), "w");
    if (!fp) err(errno, "cannot open output: %s", a.output.c_str());
    std::setvbuf(fp, nullptr, _IOFBF, 1 << 22);
    for (auto& b : bufs) std::fwrite(b.data(), 1, b.size(), fp);
    std::fclose(fp);
    std::cerr << "dist time:   " << get_sec() - t1 << " s\n";
    std::cerr << "output: " << a.output << "\n";
}

// ─────────────────────────────────────────────────────────────────────────────
//  List-mode: inverted-index pipeline (FastKMV / ProbMinHash / Kssd / SetSketch)
// ─────────────────────────────────────────────────────────────────────────────
// FastKMV indexed pipeline.  Mirrors test_FastKMV.cpp exactly:
//
//   1. Each genome's bottom-k hash list is published into a per-thread local
//      inverted index (hash -> list of genome IDs that contain it).
//   2. Local indices are merged into a single CSR-format inverted index.
//   3. computeDistances() walks each genome's keys, follows the posting list
//      for each key, and uses stamp/epoch counting to accumulate the
//      intersection size c against every other genome in O(N · K · avg_pl).
//
// Jaccard is then estimated with the standard set formula c/(s0+s1-c).  Note
// that FastKMV::jaccard() uses a slightly different estimator that caps the
// union at k (the bottom-k of A∪B), so distances from --index and no-index
// paths can differ by a small amount on the same sketch pair — both are valid
// estimators of the true Jaccard, just over different denominators.
static void run_index_fastkmv(const Args& a, const std::vector<std::string>& files) {
    const int N = static_cast<int>(files.size());
    std::vector<int>                   sketchSizes(N);
    std::vector<std::vector<uint64_t>> skKeys(N);

    const int actualThreads = std::min(a.threads, N);
    std::vector<phmap::flat_hash_map<uint64_t, std::vector<uint32_t>>> threadIdx(actualThreads);

    double t0 = get_sec();
    #pragma omp parallel num_threads(a.threads)
    {
        int tid = omp_get_thread_num();
        auto& localIdx = (tid < actualThreads) ? threadIdx[tid] : threadIdx[0];

        #pragma omp for schedule(dynamic)
        for (int t = 0; t < N; ++t) {
            gzFile fp = gzopen(files[t].c_str(), "r");
            if (!fp) continue;
            kseq_t* ks = kseq_init(fp);

            Sketch::FastKMV sk(a.fkmvK, a.kmerSize, a.seed);
            while (kseq_read(ks) >= 0)
                sk.update(ks->seq.s, ks->seq.l);

            kseq_destroy(ks);
            gzclose(fp);

            const uint64_t* regs = sk.getRegisters();
            uint32_t sz = sk.size();
            sketchSizes[t] = static_cast<int>(sz);
            skKeys[t].resize(sz);
            for (uint32_t i = 0; i < sz; ++i) {
                skKeys[t][i] = regs[i];
                localIdx[regs[i]].push_back(static_cast<uint32_t>(t));
            }
        }
    }
    double t1 = get_sec();
    std::cerr << "sketch + local index: " << t1 - t0 << " s\n";

    auto csrIdx = Sketch::buildCSRIndex<uint64_t>(threadIdx, a.threads);
    double t2 = get_sec();
    std::cerr << "build CSR index: " << t2 - t1 << " s\n";

    // Mash-distance ↔ Jaccard pruning identical to test_FastKMV.cpp:
    //   minJac = exp(-k·D) / (2 - exp(-k·D)),  minCommon = ceil(minJac · K)
    const double p_exp  = std::exp(-static_cast<double>(a.kmerSize) * a.maxDist);
    const double minJac = p_exp / (2.0 - p_exp);
    const int    minCommon = std::max(1,
        static_cast<int>(std::ceil(minJac * static_cast<double>(a.fkmvK))));
    std::cerr << "pruning: minCommon=" << minCommon << "/" << a.fkmvK
              << "  (minJac=" << minJac << ", mashD<" << a.maxDist
              << ", k=" << a.kmerSize << ")\n";

    auto setJaccard = [](int c, int s0, int s1) -> double {
        const int denom = s0 + s1 - c;
        return (denom <= 0) ? 0.0 : static_cast<double>(c) / denom;
    };
    auto minCommonFn = [minCommon](int) { return minCommon; };

    double t3 = get_sec();
    Sketch::computeDistances<uint64_t>(csrIdx, skKeys, sketchSizes, files,
        N, a.kmerSize, a.maxDist, setJaccard, minCommonFn, a.output, a.threads);
    double t4 = get_sec();
    std::cerr << "dist time: "  << t4 - t3 << " s\n";
    std::cerr << "total time: " << t4 - t0 << " s\n";
}

static void run_index_probmh(const Args& a, const std::vector<std::string>& files) {
    const int N = static_cast<int>(files.size());
    const int mSize = static_cast<int>(a.pmhM);
    std::vector<int> sketchSizes(N, mSize);
    std::vector<std::vector<uint64_t>> skKeys(N);

    const int actualThreads = std::min(a.threads, N);
    std::vector<phmap::flat_hash_map<uint64_t, std::vector<uint32_t>>> threadIdx(actualThreads);

    double t0 = get_sec();
    #pragma omp parallel num_threads(a.threads)
    {
        int tid = omp_get_thread_num();
        auto& localIdx = (tid < actualThreads) ? threadIdx[tid] : threadIdx[0];

        #pragma omp for schedule(dynamic)
        for (int t = 0; t < N; ++t) {
            Sketch::ProbMinHash4 sk(a.pmhM, a.kmerSize, a.seed, a.pmhMaxL);
            stream_seq(files[t], [&](char* seq, uint64_t len) {
                if (len < static_cast<uint64_t>(a.kmerSize)) return;
                if (a.pmhEntropy) sk.updateEntropy(seq, len);
                else              sk.update(seq, len);
            });
            sk.getInvertedIndexKeys(skKeys[t]);
            for (uint64_t key : skKeys[t])
                localIdx[key].push_back(static_cast<uint32_t>(t));
        }
    }
    double t1 = get_sec();
    std::cerr << "sketch + local index: " << t1 - t0 << " s\n";

    auto csrIdx = Sketch::buildCSRIndex<uint64_t>(threadIdx, a.threads);
    double t2 = get_sec();
    std::cerr << "build CSR index: " << t2 - t1 << " s\n";

    const double p_exp  = std::exp(-static_cast<double>(a.kmerSize) * a.maxDist);
    const double minJac = p_exp / (2.0 - p_exp);
    const int minCommon = std::max(1, static_cast<int>(std::ceil(minJac * a.pmhM)));
    std::cerr << "pruning: minCommon=" << minCommon << "/" << a.pmhM
              << "  (minJac=" << minJac << ", mashD<" << a.maxDist
              << ", k=" << a.kmerSize << ")\n";

    auto jaccardFn = [mSize](int c, int, int) -> double {
        return Sketch::ProbMinHash4::jaccardFromCommon(c, static_cast<uint32_t>(mSize));
    };
    auto minCommonFn = [minCommon](int) { return minCommon; };

    double t3 = get_sec();
    Sketch::computeDistances<uint64_t>(csrIdx, skKeys, sketchSizes, files,
        N, a.kmerSize, a.maxDist, jaccardFn, minCommonFn, a.output, a.threads);
    double t4 = get_sec();
    std::cerr << "dist time: "  << t4 - t3 << " s\n";
    std::cerr << "total time: " << t4 - t0 << " s\n";
}

static void run_index_kssd(const Args& a, const std::vector<std::string>& files) {
    const int N = static_cast<int>(files.size());
    Sketch::kssd_parameter_t P(/*half_k*/ a.kmerSize / 2, /*half_subk*/ 6, a.kssdDrlevel);
    const int kmer = static_cast<int>(P.kmer_size);
    const bool use64 = (P.half_k - P.drlevel) > 8;

    std::vector<std::vector<uint64_t>> skKeys(N);
    std::vector<int> sketchSizes(N, 0);

    const int actualThreads = std::min(a.threads, N);
    std::vector<phmap::flat_hash_map<uint64_t, std::vector<uint32_t>>> threadIdx(actualThreads);

    double t0 = get_sec();
    #pragma omp parallel num_threads(a.threads)
    {
        int tid = omp_get_thread_num();
        auto& localIdx = (tid < actualThreads) ? threadIdx[tid] : threadIdx[0];

        #pragma omp for schedule(dynamic)
        for (int t = 0; t < N; ++t) {
            Sketch::Kssd sk(P);
            stream_seq(files[t], [&](char* seq, uint64_t /*len*/) { sk.update(seq); });

            if (use64) {
                auto vals = sk.storeHashes64();
                std::sort(vals.begin(), vals.end());
                vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
                skKeys[t] = std::move(vals);
            } else {
                auto vals32 = sk.storeHashes();
                std::sort(vals32.begin(), vals32.end());
                vals32.erase(std::unique(vals32.begin(), vals32.end()), vals32.end());
                skKeys[t].reserve(vals32.size());
                for (uint32_t h : vals32)
                    skKeys[t].push_back(static_cast<uint64_t>(h));
            }
            sketchSizes[t] = static_cast<int>(skKeys[t].size());
            for (uint64_t h : skKeys[t])
                localIdx[h].push_back(static_cast<uint32_t>(t));
        }
    }
    double t1 = get_sec();
    std::cerr << "sketch + local index: " << t1 - t0 << " s\n";

    auto csrIdx = Sketch::buildCSRIndex<uint64_t>(threadIdx, a.threads);
    double t2 = get_sec();
    std::cerr << "build CSR index: " << t2 - t1 << " s\n";

    const double p_exp  = std::exp(-static_cast<double>(kmer) * a.maxDist);
    const double minJac = p_exp / (2.0 - p_exp);
    std::cerr << "pruning: per-sketch minCommon=ceil(" << minJac << " * |s|)"
              << "  (mashD<" << a.maxDist << ", k=" << kmer << ")\n";

    auto jaccardFn = [](int c, int s0, int s1) -> double {
        const int denom = s0 + s1 - c;
        return (denom <= 0) ? 0.0 : static_cast<double>(c) / denom;
    };
    auto minCommonFn = [minJac](int s0) {
        return std::max(1, static_cast<int>(std::ceil(minJac * s0)));
    };

    double t3 = get_sec();
    Sketch::computeDistances<uint64_t>(csrIdx, skKeys, sketchSizes, files,
        N, kmer, a.maxDist, jaccardFn, minCommonFn, a.output, a.threads);
    double t4 = get_sec();
    std::cerr << "dist time: "  << t4 - t3 << " s\n";
    std::cerr << "total time: " << t4 - t0 << " s\n";
}

static void run_index_setsketch(const Args& a, const std::vector<std::string>& files) {
    const int N = static_cast<int>(files.size());
    Sketch::SetSketch proto(a.ssBits, a.ssBase, a.ssA);
    const int    m      = proto.getM();
    const double factor = proto.getFactor();
    double bip_buf[64];
    std::memcpy(bip_buf, proto.getBaseInvPow(), 64 * sizeof(double));
    const double* bip = bip_buf;

    static const int WITNESS_STRIDE = 4;
    const int witnessesPerSketch = m / WITNESS_STRIDE;

    std::vector<double>  sizes(N, 0.0);
    std::vector<uint8_t> flat_cores(static_cast<size_t>(N) * m, 0);
    std::vector<std::vector<uint64_t>> skKeys(N);

    double t0 = get_sec();
    #pragma omp parallel for num_threads(a.threads) schedule(dynamic)
    for (int t = 0; t < N; ++t) {
        Sketch::SetSketch sk(a.ssBits, a.ssBase, a.ssA);
        stream_seq(files[t], [&](char* seq, uint64_t len) {
            sk.update(seq, static_cast<size_t>(len));
        });
        sizes[t] = sk.cardinality();
        std::memcpy(&flat_cores[static_cast<size_t>(t) * m],
                    sk.getCore().data(), m);
        const uint64_t* wit = sk.getWitnesses().data();
        skKeys[t].reserve(witnessesPerSketch);
        for (int p = 0; p < m; p += WITNESS_STRIDE)
            if (wit[p] != 0) skKeys[t].push_back(wit[p]);
    }
    double t1 = get_sec();
    std::cerr << "sketch time: " << t1 - t0 << " s\n";

    static const int TAIL_STEP = 64;
    const int nCP = m / TAIL_STEP + 1;
    std::vector<double> tailSums(static_cast<size_t>(N) * nCP, 0.0);

    #pragma omp parallel for num_threads(a.threads) schedule(static)
    for (int t = 0; t < N; ++t) {
        const uint8_t* core = &flat_cores[static_cast<size_t>(t) * m];
        double* ts = &tailSums[static_cast<size_t>(t) * nCP];
        ts[nCP - 1] = 0.0;
        for (int cp = nCP - 2; cp >= 0; --cp) {
            const int start = cp * TAIL_STEP;
            const int end   = std::min(start + TAIL_STEP, m);
            double blockSum = 0.0;
            for (int kk = start; kk < end; ++kk) blockSum += bip[core[kk]];
            ts[cp] = ts[cp + 1] + blockSum;
        }
    }

    const int actualThreads = std::min(a.threads, N);
    std::vector<phmap::flat_hash_map<uint64_t, std::vector<uint32_t>>> threadIdx(actualThreads);

    #pragma omp parallel num_threads(a.threads)
    {
        int tid = omp_get_thread_num();
        if (tid < actualThreads) {
            auto& localIdx = threadIdx[tid];
            #pragma omp for schedule(static)
            for (int t = 0; t < N; ++t)
                for (uint64_t key : skKeys[t])
                    localIdx[key].push_back(static_cast<uint32_t>(t));
        }
    }

    auto csrIdx = Sketch::buildCSRIndex<uint64_t>(threadIdx, a.threads);
    double t2 = get_sec();
    std::cerr << "build CSR index + tailSums: " << t2 - t1 << " s\n";

    const double p_exp  = std::exp(-static_cast<double>(a.kmerSize) * a.maxDist);
    const double minJac = p_exp / (2.0 - p_exp);
    const double expectedOverlap = static_cast<double>(witnessesPerSketch) * minJac;
    const double sd = std::sqrt(std::max(0.0, expectedOverlap * (1.0 - minJac)));
    const int minCommon = std::max(1, static_cast<int>(std::floor(expectedOverlap - 8.0 * sd)));
    std::cerr << "pruning: minCommon=" << minCommon << "/" << witnessesPerSketch
              << "  (minJac=" << minJac << ", mashD<" << a.maxDist
              << ", k=" << a.kmerSize << ")\n";

    auto exactJaccardFn = [&](int i, int j) -> double {
        const double si = sizes[i], sj = sizes[j];
        if (si > 0.0 && sj / si < minJac) return -1.0;
        return Sketch::SetSketch::jaccardFromCoresBatch(
            &flat_cores[static_cast<size_t>(i) * m],
            &flat_cores[static_cast<size_t>(j) * m],
            m, bip, factor, si, sj, minJac,
            &tailSums[static_cast<size_t>(i) * nCP],
            &tailSums[static_cast<size_t>(j) * nCP],
            TAIL_STEP);
    };
    auto minCommonFn = [minCommon](int) { return minCommon; };

    double t3 = get_sec();
    Sketch::computeDistancesExact<uint64_t>(csrIdx, skKeys, files,
        N, a.kmerSize, a.maxDist, exactJaccardFn, minCommonFn, a.output, a.threads);
    double t4 = get_sec();
    std::cerr << "dist time: "  << t4 - t3 << " s\n";
    std::cerr << "total time: " << t4 - t0 << " s\n";
}

// ─────────────────────────────────────────────────────────────────────────────
//  Main
// ─────────────────────────────────────────────────────────────────────────────
int main(int argc, char* argv[]) {
    Args a;
    if (parse_args(argc, argv, a) != 0) return 1;

    if (!a.listMode) return run_pair(a);

    auto files = load_file_list(a.inputs[0]);
    if (files.empty()) { std::cerr << "ERROR: empty file list\n"; return 1; }
    std::cerr << "===== total files: " << files.size() << "  algo=";
    switch (a.algo) {
        case Algo::MINHASH:     std::cerr << "MinHash";     break;
        case Algo::KSSD:        std::cerr << "Kssd";        break;
        case Algo::HLL:         std::cerr << "HLL";         break;
        case Algo::PROBMINHASH: std::cerr << "ProbMinHash"; break;
        case Algo::BINDASH:     std::cerr << "BinDash";     break;
        case Algo::SETSKETCH:   std::cerr << "SetSketch";   break;
        case Algo::FASTKMV:     std::cerr << "FastKMV";     break;
        default: break;
    }
    std::cerr << (a.useIndex ? "  (inverted index)\n" : "  (all-pairs)\n");

    if (a.useIndex) {
        switch (a.algo) {
            case Algo::FASTKMV:     run_index_fastkmv(a, files);    return 0;
            case Algo::PROBMINHASH: run_index_probmh(a, files);     return 0;
            case Algo::KSSD:        run_index_kssd(a, files);       return 0;
            case Algo::SETSKETCH:   run_index_setsketch(a, files);  return 0;
            default:
                std::cerr << "WARNING: --index not supported for this algorithm; "
                             "falling back to all-pairs\n";
                break;
        }
    }

    // All-pairs (no index) for every algorithm.
    const int kmer = a.kmerSize;
    switch (a.algo) {
    case Algo::MINHASH:
        list_allpairs<Sketch::MinHash>(a, files, build_minhash,
            [](Sketch::MinHash* x, Sketch::MinHash* y) { return x->distance(y); });
        break;

    case Algo::HLL:
        list_allpairs<Sketch::HyperLogLog>(a, files, build_hll,
            [kmer](Sketch::HyperLogLog* x, Sketch::HyperLogLog* y) {
                return mash_distance(x->jaccard_index(*y), kmer);
            });
        break;

    case Algo::KSSD:
        list_allpairs<Sketch::Kssd>(a, files, build_kssd,
            [](Sketch::Kssd* x, Sketch::Kssd* y) { return x->distance(y); });
        break;

    case Algo::PROBMINHASH:
        list_allpairs<Sketch::ProbMinHash4>(a, files, build_probminhash,
            [kmer](Sketch::ProbMinHash4* x, Sketch::ProbMinHash4* y) {
                return mash_distance(x->jaccard(*y), kmer);
            });
        break;

    case Algo::BINDASH:
        list_allpairs<Sketch::BinDash>(a, files, build_bindash,
            [kmer](Sketch::BinDash* x, Sketch::BinDash* y) {
                return mash_distance(x->jaccard(*y), kmer);
            });
        break;

    case Algo::SETSKETCH:
        list_allpairs<Sketch::SetSketch>(a, files, build_setsketch,
            [kmer](Sketch::SetSketch* x, Sketch::SetSketch* y) {
                return mash_distance(x->jaccard_index(*y), kmer);
            });
        break;

    case Algo::FASTKMV:
        list_allpairs<Sketch::FastKMV>(a, files, build_fastkmv,
            [kmer](Sketch::FastKMV* x, Sketch::FastKMV* y) {
                return mash_distance(x->jaccard(*y), kmer);
            });
        break;

    default: return 1;
    }

    return 0;
}
