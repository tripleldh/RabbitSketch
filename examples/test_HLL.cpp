#include "Sketch.h"
//#include <iostream>
#include <sys/time.h>
#include <zlib.h>
#include "kseq.h"
#include <vector>
#include <math.h>
#include <random>
#include <fstream>
#include <err.h>
#include <sys/stat.h>
#include <omp.h>
#include <sstream>
#include <atomic>
#include <algorithm>
#include "common.h"
#if defined(__GNUC__) && defined(_OPENMP)
# include <parallel/algorithm>
#endif
using namespace std;

typedef struct fileInfo
{
  string fileName;
} fileInfo_t;

KSEQ_INIT(gzFile, gzread)

//  double get_sec(){
//    struct timeval tv;
//    gettimeofday(&tv, NULL);
//    return tv.tv_sec + (double)tv.tv_usec/1000000;
//  }


int main(int argc, char* argv[])
{
	if(argc < 4){
		cerr << "usage: " << argv[0] << " <file_list> <dist_threshold> <threads>" << endl;
		return 1;
	}
  string inputFile = argv[1];
  double thres = stod(argv[2]);
  int numThreads = stoi(argv[3]);
  ifstream fs(inputFile);
  if(!fs){
    err(errno, "cannot open the inputFile: %s\n", inputFile.c_str());
  }
  vector<fileInfo_t> fileList;
  uint64_t totalSize = 0;
  string fileName;
  while(getline(fs, fileName)){
    struct stat cur_stat;
    stat(fileName.c_str(), &cur_stat);
    fileInfo_t tmpF;
    tmpF.fileName = fileName;
    fileList.push_back(tmpF);
  }
  double t1 = get_sec();
  vector<string> fileArr;
  for(size_t i = 0; i < fileList.size(); i++){
    fileArr.push_back(fileList[i].fileName);
  }
  int small_file_number = fileArr.size();
  vector<Sketch::HyperLogLog> vhlog;
  cerr << "=====total small files: " << small_file_number << endl;
	vector<string> resFileName;
	#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
  for(size_t t = 0; t <small_file_number; t++)
  {
    gzFile fp1;
    kseq_t * ks1;
    fp1 = gzopen(fileArr[t].c_str(), "r");
    if(fp1 == NULL){
      err(errno, "cannot open the genome file: %s\n", fileArr[t].c_str());
    }
    ks1 = kseq_init(fp1);
    static const size_t BITS = 10; //24
    Sketch::HyperLogLog hll(BITS);
		while(1){
      int length = kseq_read(ks1);
      if(length < 0){
        break;
 			}

     	hll.update(ks1->seq.s);
    }
			#pragma omp critical
      {
        vhlog.push_back(hll);
				resFileName.push_back(fileArr[t]);
      }
    gzclose(fp1);
    kseq_destroy(ks1);
  }
  double t2 = get_sec();
  cerr << "sketch time is: " << t2 - t1 << endl;

	vector<string> thread_bufs(numThreads);

	const int n = (int)vhlog.size();
	cerr << "vhlog size is: " << n << endl;

	cerr << "LSH index breakdown (single-thread steps marked):" << endl;

	// ── Phase 1: pre-compute cardinalities (parallel) ────────────────────────
	vector<double> sizes(n);
	#pragma omp parallel for num_threads(numThreads) schedule(static)
	for(int i = 0; i < n; i++)
		sizes[i] = vhlog[i].cardinality();
	double t_phase1 = get_sec();
	cerr << "  Phase 1 cardinality (parallel): " << t_phase1 - t2 << " s" << endl;

	// ── Phase 2: LSH banding ─────────────────────────────────────────────────
	const int m     = (n > 0) ? (int)vhlog[0].getCore().size() : 0;
	const int BANDS = 128;
	const int ROWS  = (m > 0) ? m / BANDS : 8;

	auto band_hash = [](const uint8_t* data, int len) -> uint32_t {
		uint32_t h = 2166136261u;
		for(int i = 0; i < len; i++){ h ^= data[i]; h *= 16777619u; }
		return h;
	};

	vector<pair<uint64_t,int>> band_entries((size_t)n * BANDS);
	#pragma omp parallel for num_threads(numThreads) schedule(static)
	for(int i = 0; i < n; i++){
		const auto& core = vhlog[i].getCore();
		for(int b = 0; b < BANDS; b++){
			uint32_t h = band_hash(core.data() + b * ROWS, ROWS);
			band_entries[(size_t)i * BANDS + b] = { ((uint64_t)(uint32_t)b << 32) | (uint32_t)h, i };
		}
	}
	double t_build = get_sec();
	cerr << "  Phase 2 build band_entries (parallel): " << t_build - t_phase1 << " s" << endl;

	sort(band_entries.begin(), band_entries.end());
	double t_sort_bands = get_sec();
	cerr << "  Phase 2 sort band_entries [SINGLE]: " << t_sort_bands - t_build << " s" << endl;

	// Build key-group boundaries (one pass), then emit pairs in parallel per group.
	vector<size_t> group_start;
	group_start.push_back(0);
	for(size_t i = 1; i < band_entries.size(); i++)
		if(band_entries[i].first != band_entries[i-1].first)
			group_start.push_back(i);
	group_start.push_back(band_entries.size());
	const int n_bands = (int)group_start.size() - 1;

	vector<vector<pair<int,int>>> thread_cands((size_t)numThreads);

	#pragma omp parallel num_threads(numThreads)
	{
		int tid = omp_get_thread_num();
		auto& local = thread_cands[tid];
		local.clear();
		#pragma omp for schedule(dynamic)
		for(int g = 0; g < n_bands; g++){
			size_t s = group_start[g], e = group_start[g+1];
			for(size_t a = s; a < e; a++)
				for(size_t b = a+1; b < e; b++){
					int ia = band_entries[a].second, ib = band_entries[b].second;
					local.emplace_back(min(ia, ib), max(ia, ib));
				}
		}
	}
	band_entries.clear();

	vector<size_t> prefix((size_t)numThreads + 1);
	prefix[0] = 0;
	for(int t = 0; t < numThreads; t++) prefix[t+1] = prefix[t] + thread_cands[t].size();
	size_t total_cands = prefix[numThreads];
	vector<pair<int,int>> candidates(total_cands);
	#pragma omp parallel for num_threads(numThreads)
	for(int t = 0; t < numThreads; t++)
		copy(thread_cands[t].begin(), thread_cands[t].end(), candidates.begin() + prefix[t]);
	thread_cands.clear();
	thread_cands.shrink_to_fit();
	double t_scan = get_sec();
	cerr << "  Phase 2 scan→candidates (parallel): " << t_scan - t_sort_bands << " s" << endl;

#if defined(__GNUC__) && defined(_OPENMP)
	__gnu_parallel::sort(candidates.begin(), candidates.end());
#else
	sort(candidates.begin(), candidates.end());
#endif
	candidates.erase(unique(candidates.begin(), candidates.end()), candidates.end());
	double t_lsh = get_sec();
	cerr << "  Phase 2 sort+dedup candidates: " << t_lsh - t_scan << " s" << endl;

	const long long total_pairs = (long long)n * (n - 1) / 2;
	cerr << "LSH candidates: " << candidates.size()
	     << " / " << total_pairs << " total pairs"
	     << "  (reduction: " << 100.0*(1.0-(double)candidates.size()/total_pairs) << "%)" << endl;
	cerr << "LSH index time (total): " << t_lsh - t2 << " s" << endl;

	// ── Phase 3: verify candidates in parallel ────────────────────────────────
	// dist = Jaccard distance (1 - J); thres in [0,1]: 0 = identical, 1 = disjoint
	const double min_jaccard = 1.0 - thres;

	std::atomic<long long> cnt_size_filtered{0};
	std::atomic<long long> cnt_exact{0};

	#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
	for(int c = 0; c < (int)candidates.size(); c++){
		int i = candidates[c].first, j = candidates[c].second;
		int tid = omp_get_thread_num();

		// Guard 1: size ratio  (Jaccard ≤ min/max cardinality)
		// HLL cardinality has ~3.1% relative error (1/sqrt(1024)); two estimates
		// can differ by up to ~6% even for identical sets.  Apply a 7% safety
		// margin to avoid false negatives near the threshold.
		double si = sizes[i], sj = sizes[j];
		if(si <= 0 || sj <= 0) { cnt_size_filtered++; continue; }
		if(min(si, sj) / max(si, sj) < min_jaccard * 0.93) { cnt_size_filtered++; continue; }

		// Guard 2: full Ertl joint MLE
		cnt_exact++;
		double dist = vhlog[i].distance(vhlog[j]);
		if(dist < thres) {
			char line[4096];
			int len = snprintf(line, sizeof(line), "%s\t%s\t%lf\n",
			                   resFileName[i].c_str(), resFileName[j].c_str(), dist);
			thread_bufs[tid].append(line, len);
		}
	}

	// Flush all buffers to a single output file
	system("mkdir -p res_dir");
	FILE* fp_out = fopen("res_dir/res.dist.HLL", "w");
	for(int t = 0; t < numThreads; t++)
		fwrite(thread_bufs[t].data(), 1, thread_bufs[t].size(), fp_out);
	fclose(fp_out);

	double t3 = get_sec();
	cerr << "dist time is: " << t3 - t_lsh << " s" << endl;
	cerr << "  candidates:    " << (long long)candidates.size() << endl;
	cerr << "  size-filtered: " << cnt_size_filtered.load() << endl;
	cerr << "  exact computed:" << cnt_exact.load() << endl;

}




















