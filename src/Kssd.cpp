#include <omp.h>
#include "Sketch.h"
#include <cstdio>
#include <math.h>
#include <sys/stat.h>
#include <immintrin.h>
#include <cstring>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <atomic>
#include <stdint.h>
#include <queue>
#include <algorithm>
#include "robin_hood.h"
//#include "Kssd.h"
#include <err.h>
#include <fstream>
#include <cstdlib>
#include <unordered_map>
#include <unordered_set>
#include "common.h"
#define COMPONENT_SZ 7
#define MIN_SUBCTX_DIM_SMP_SZ 4096
#define _64MASK 0xffffffffffffffffLLU
#define CTX_SPC_USE_L 8
#define LD_FCTR 0.6

using namespace std;

#define H1(K, HASH_SZ) ((K) % (HASH_SZ))
#define H2(K, HASH_SZ) (1 + (K) % ((HASH_SZ)-1))
#define HASH(K, I, HASH_SZ) ((H1(K, HASH_SZ) + I * H2(K, HASH_SZ)) % HASH_SZ)

#define PATHLEN 256

using Sketch::sketchInfo_t;
bool isSketchFile(string inputFile){
	auto startPos = inputFile.find_last_of('.');
	if(startPos == string::npos) return false;
	string suffixName = inputFile.substr(startPos+1);
	if(suffixName == "sketch")	return true;
	else return false;
}

namespace Sketch{
void saveSketches(vector<Sketch::sketch_t>& sketches, Sketch::sketchInfo_t& info, string outputFile){
	bool use64 = info.half_k - info.drlevel > 8 ? true : false;
	FILE * fp = fopen(outputFile.c_str(), "w+");
	int sketchNumber = sketches.size();
	info.genomeNumber = sketchNumber;
	info.id = (info.half_k << 8) + (info.half_subk << 4) + info.drlevel;
	fwrite(&info, sizeof(Sketch::sketchInfo_t), 1, fp);

	int * genomeNameSize = new int[sketchNumber];
	int * hashSetSize = new int[sketchNumber];
	//uint64_t totalNumber = 0;
	//uint64_t totalLength = 0;
	for(int i = 0; i < sketchNumber; i++){
		genomeNameSize[i] = sketches[i].fileName.length();
		if(use64)
			hashSetSize[i] = sketches[i].hashSet64.size();
		else
			hashSetSize[i] = sketches[i].hashSet.size();
		//totalNumber += hashSetSize[i];
		//totalLength += genomeNameSize[i];
	}
	//cerr << "the sketches size is: " << sketchNumber << endl;
	//cerr << "the total hash number is: " << totalNumber << endl;
	//cerr << "the total name length is: " << totalLength << endl;
	//fwrite(&parameter, sizeof(kssd_parameter_t), 1, fp);
	//fwrite(&sketchNumber, sizeof(int), 1, fp);
	fwrite(genomeNameSize, sizeof(int), sketchNumber, fp);
	fwrite(hashSetSize, sizeof(int), sketchNumber, fp);
	for(int i = 0; i < sketchNumber; i++){
		const char * namePoint = sketches[i].fileName.c_str();
		fwrite(namePoint, sizeof(char), genomeNameSize[i], fp);
		if(use64){
			uint64_t * curPoint = sketches[i].hashSet64.data();
			fwrite(curPoint, sizeof(uint64_t), hashSetSize[i], fp);
		}
		else{
			uint32_t * curPoint = sketches[i].hashSet.data();
			fwrite(curPoint, sizeof(uint32_t), hashSetSize[i], fp);
		}
	}
	fclose(fp);
	delete [] genomeNameSize;
	delete [] hashSetSize;

}

void readSketches(vector<Sketch::sketch_t>& sketches, Sketch::sketchInfo_t& info, string inputFile){
	FILE * fp = fopen(inputFile.c_str(), "rb");
	if(!fp){
		fprintf(stderr, "ERROR: readSketches(), cannot open the file: %s\n", inputFile.c_str());
		exit(1);
	}
	int read_sketch_info = fread(&info, sizeof(Sketch::sketchInfo_t), 1, fp);
	bool use64 = info.half_k - info.drlevel > 8 ? true : false;
	int sketchNumber = info.genomeNumber;
	//fread(&sketchNumber, sizeof(int), 1, fp);
	//cerr << "sketchNumber is: " << sketchNumber << endl;
	int * genomeNameSize = new int[sketchNumber];
	int * hashSetSize = new int[sketchNumber];
	int read_genome_name_size = fread(genomeNameSize, sizeof(int), sketchNumber, fp);
	int read_hash_set_size = fread(hashSetSize, sizeof(int), sketchNumber, fp);
	if(read_sketch_info != 1 || read_genome_name_size != sketchNumber || read_hash_set_size != sketchNumber){
		cerr << "ERROR: readSketches(), mismatched read sketch_info, genome_name_sizse, hash_set_size" << endl;
		exit(1);
	}

	//uint64_t totalNumber = 0;
	//uint64_t totalLength = 0;

	int maxNameLength = 1000;
	char * curName = new char[maxNameLength+1];
	int maxHashSize = 1 << 24;
	uint32_t * curPoint = new uint32_t[maxHashSize];
	uint64_t * curPoint64 = new uint64_t[maxHashSize];
	for(int i = 0; i < sketchNumber; i++){
		//read the genome name.
		int curLength = genomeNameSize[i];
		if(curLength > maxNameLength){
			maxNameLength = curLength;
			curName = new char[maxNameLength+1];
		}
		//totalLength += curLength;
		int nameLength = fread(curName, sizeof(char), curLength, fp);
		if(nameLength != curLength){
			cerr << "error: the read nameLength is not equal to the saved nameLength, exit!" << endl;
			exit(0);
		}
		string genomeName;
		genomeName.assign(curName, curName + curLength);

		//read the hash values in each sketch
		int curSize = hashSetSize[i];
		//totalNumber += curSize;
		if(curSize > maxHashSize){
			maxHashSize = curSize;
			if(use64)
				curPoint64 = new uint64_t[maxHashSize];
			else
				curPoint = new uint32_t[maxHashSize];
		}

		Sketch::sketch_t s;
		s.fileName = genomeName;
		s.id = i;
		if(use64){
			int hashSize = fread(curPoint64, sizeof(uint64_t), curSize, fp);
			if(hashSize != curSize){
				cerr << "ERROR: readSketches(), the read hashNumber for a sketch is not equal to the saved hashNumber, exit" << endl;
				exit(1);
			}
			vector<uint64_t> curHashSet64(curPoint64, curPoint64 + curSize);
			s.hashSet64=curHashSet64;
		}
		else{
			int hashSize = fread(curPoint, sizeof(uint32_t), curSize, fp);
			if(hashSize != curSize){
				cerr << "ERROR: readSketches(), the read hashNumber for a sketch is not equal to the saved hashNumber, exit" << endl;
				exit(0);
			}
			vector<uint32_t> curHashSet(curPoint, curPoint + curSize);
			s.hashSet=curHashSet;
		}
		sketches.push_back(s);
	}
	delete [] curPoint;
	delete [] curName;
	fclose(fp);
	//cerr << "the total number of hash value of: " << inputFile << " is: " << totalNumber << endl;
	//cerr << "the total length of genome name of: " << inputFile << " is: " << totalLength << endl;

}



void transSketches(vector<Sketch::sketch_t>& sketches, Sketch::sketchInfo_t& info, string dictFile, string indexFile, int numThreads){
	double t0 = get_sec();
	int half_k = info.half_k;
	int drlevel = info.drlevel;
	bool use64 = half_k - drlevel > 8 ? true : false;
	if(use64)
		cerr << "transSketches: use64" << endl;
	else
		cerr << "transSketches: not use64" << endl;

	if(use64){
		double t0 = get_sec();
		//size_t dict_size = (1LLU << (4*(half_k-drlevel))) / 64;
		//uint64_t* dict = (uint64_t*)malloc(dict_size * sizeof(uint64_t));
		//memset(dict, 0, dict_size * sizeof(uint64_t));
		robin_hood::unordered_map<uint64_t, vector<uint32_t>> hash_map_arr;
		//std::unordered_map<uint64_t, vector<uint32_t>> hash_map_arr;
		//std::map<uint64_t, vector<uint32_t>> hash_map_arr;
		for(size_t i = 0; i < sketches.size(); i++){
			//#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
			for(size_t j = 0; j < sketches[i].hashSet64.size(); j++){
				uint64_t cur_hash = sketches[i].hashSet64[j];
				//cerr << cur_hash << endl;
				//dict[cur_hash/64] |= (0x8000000000000000LLU >> (cur_hash % 64));
				hash_map_arr.insert({cur_hash, vector<uint32_t>()});
				hash_map_arr[cur_hash].push_back(i);
			}
		}
		//string filterFile = indexFile + ".filter";
		//FILE* fp_filter = fopen(filterFile.c_str(), "w+");
		//fwrite(dict, sizeof(uint64_t), dict_size, fp_filter);
		//fclose(fp_filter);
		//free(dict);
		double t1 = get_sec();
#ifdef Timer_inner
		cerr << "the time of generate the bloom dictionary and hash_map_arr is: " << t1 - t0 << endl;
#endif
		size_t hash_number = hash_map_arr.size();
		cerr << "the hash_number is: " << hash_number << endl;
		size_t total_size = 0;
		uint64_t* hash_arr = (uint64_t*)malloc(hash_number * sizeof(uint64_t));
		uint32_t* hash_size_arr = (uint32_t*)malloc(hash_number * sizeof(uint32_t));
		FILE* fp_dict = fopen(dictFile.c_str(), "w+");
		if(!fp_dict){
			cerr << "ERROR: transSketches, cannot open dictFile: " << dictFile << endl;
			exit(1);
		}
		size_t cur_id = 0;
		for(auto x : hash_map_arr){
			hash_arr[cur_id] = x.first;
			fwrite(x.second.data(), sizeof(uint32_t), x.second.size(), fp_dict);
			hash_size_arr[cur_id] = x.second.size();
			total_size += x.second.size();
			cur_id++;
		}
		cerr << "the total size is: " << total_size << endl;
		fclose(fp_dict);
		double t2 = get_sec();
#ifdef Timer_inner
		cerr << "the time of writing dictFile is: " << t2 - t1 << endl;
#endif

		FILE* fp_index = fopen(indexFile.c_str(), "w+");
		if(!fp_index){
			cerr << "ERROR: transSketches, cannot open indexFile: " << indexFile << endl;
			exit(1);
		}
		fwrite(&hash_number, sizeof(size_t), 1, fp_index);
		fwrite(hash_arr, sizeof(uint64_t), hash_number, fp_index);
		fwrite(hash_size_arr, sizeof(uint32_t), hash_number, fp_index);
		fclose(fp_index);
		double t3 = get_sec();
#ifdef Timer_inner
		cerr << "the time of writing indexFile is: " << t3 - t2 << endl;
#endif
	}
	else{
		size_t hashSize = 1LLU << (4 * (half_k - drlevel));
		vector<vector<uint32_t>> hashMapId;
		for(size_t i = 0; i < hashSize; i++){
			hashMapId.push_back(vector<uint32_t>());
		}
		uint32_t* offsetArr = (uint32_t*)calloc(hashSize, sizeof(uint32_t));

		cerr << "the hashSize is: " << hashSize << endl;
		for(size_t i = 0; i < sketches.size(); i++){
#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
			for(size_t j = 0; j < sketches[i].hashSet.size(); j++){
				uint32_t hash = sketches[i].hashSet[j];
				hashMapId[hash].push_back(i);
			}
		}
		double tt0 = get_sec();
#ifdef Timer_inner
		cerr << "the time of generate the idx by multiple threads are: " << tt0 - t0 << endl;
#endif

		FILE * fp0 = fopen(dictFile.c_str(), "w+");
		uint64_t totalIndex = 0;
		for(size_t hash = 0; hash < hashSize; hash++){
			offsetArr[hash] = 0;
			if(hashMapId[hash].size() != 0){
				fwrite(hashMapId[hash].data(), sizeof(uint32_t), hashMapId[hash].size(), fp0);
				totalIndex += hashMapId[hash].size();
				offsetArr[hash] = hashMapId[hash].size();
			}
		}
		fclose(fp0);

		double t1 = get_sec();
#ifdef Timer_inner
		cerr << "the time of merge multiple idx into final hashMap is: " << t1 - tt0 << endl;
#endif

		FILE * fp1 = fopen(indexFile.c_str(), "w+");
		fwrite(&hashSize, sizeof(size_t), 1, fp1);
		fwrite(&totalIndex, sizeof(uint64_t), 1, fp1);
		fwrite(offsetArr, sizeof(uint32_t), hashSize, fp1);
		double t2 = get_sec();
		fclose(fp1);
#ifdef Timer_inner
		cerr << "the time of write output file is: " << t2 - t1 << endl;
#endif
	}

	//cerr << "the hashSize is: " << hashSize << endl;
	//cerr << "the totalIndex is: " << totalIndex << endl;
}


void printInfos(vector<Sketch::sketch_t>& sketches, string outputFile){
	FILE * fp = fopen(outputFile.c_str(), "w+");
	fprintf(fp, "the number of sketches are: %lu\n", sketches.size());
	for(size_t i = 0; i < sketches.size(); i++){
		fprintf(fp, "%s\t%lu\n", sketches[i].fileName.c_str(), sketches[i].hashSet.size());
	}
	fclose(fp);
}

void printSketches(vector<Sketch::sketch_t>& sketches, string outputFile){
	FILE * fp = fopen(outputFile.c_str(), "w+");
	fprintf(fp, "the number of sketches are: %lu\n", sketches.size());
	for(size_t i = 0; i < sketches.size(); i++){
		fprintf(fp, "%s\t%lu\n", sketches[i].fileName.c_str(), sketches[i].hashSet.size());
		for(size_t j = 0; j < sketches[i].hashSet.size(); j++){
			fprintf(fp, "%u\t", sketches[i].hashSet[j]);
			if(j % 10 == 9) fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

//bool sketchFastaFile(string inputFile, bool isQuery, int numThreads, kssd_parameter_t parameter, vector<sketch_t>& sketches, sketchInfo_t& info, string outputFile) {
//    //cerr << "run the sketchFastaFile " << endl;
//    int half_k = parameter.half_k;
//    int half_subk = parameter.half_subk;
//    int drlevel = parameter.drlevel;
//    int rev_add_move = parameter.rev_add_move;
//    int half_outctx_len = parameter.half_outctx_len;
//    int* shuffled_dim = parameter.shuffled_dim;
//    int dim_start = parameter.dim_start;
//    int dim_end = parameter.dim_end;
//    int kmer_size = parameter.kmer_size;
//    uint64_t tupmask = parameter.tupmask;
//    uint64_t undomask0 = parameter.undomask0;
//    uint64_t undomask1 = parameter.undomask1;
//    uint64_t domask = parameter.domask;
//
//    bool use64 = half_k - drlevel > 8 ? true : false;
//
//    int dim_size = 1 << 4 * (half_k - half_outctx_len);
//    robin_hood::unordered_map<uint32_t, int> shuffled_map;
//    for (int t = 0; t < dim_size; t++) {
//        if (shuffled_dim[t] < dim_end && shuffled_dim[t] >= dim_start) {
//            shuffled_map.insert({t, shuffled_dim[t]});
//        }
//    }
//
//    ifstream fs(inputFile);
//    if (!fs) {
//        err(errno, "cannot open the inputFile: %s\n", inputFile.c_str());
//    }
//    vector<fileInfo_t> fileList;
//    uint64_t totalSize = 0;
//    string fileName;
//    while (getline(fs, fileName)) {
//        struct stat cur_stat;
//        stat(fileName.c_str(), &cur_stat);
//        uint64_t curSize = cur_stat.st_size;
//        totalSize += curSize;
//        fileInfo_t tmpF;
//        tmpF.fileName = fileName;
//        tmpF.fileSize = curSize;
//        fileList.push_back(tmpF);
//    }
//    std::sort(fileList.begin(), fileList.end(), cmpFile);
//
//    // All files are considered small files now, so no need for big/small file differentiation
//    int small_file_number = fileList.size();
//    int process_bar_size = get_progress_bar_size(small_file_number);
//    cerr << "=====total small files: " << small_file_number << endl;
//
//    #pragma omp parallel for num_threads(numThreads) schedule(dynamic)
//    for (size_t t = 0; t < small_file_number; t++) {
//        //int tid = omp_get_thread_num();
//        sketch_t tmpSketch;
//        gzFile fp1;
//        kseq_t* ks1;
//        fp1 = gzopen(fileList[t].fileName.c_str(), "r");
//        if (fp1 == NULL) {
//            err(errno, "cannot open the genome file: %s\n", fileList[t].fileName.c_str());
//        }
//        ks1 = kseq_init(fp1);
//        uint64_t totalLength = 0;
//        tmpSketch.fileName = fileList[t].fileName;
//
//        unordered_set<uint32_t> hashValueSet;
//        unordered_set<uint64_t> hashValueSet64;
//        while (1) {
//            int length = kseq_read(ks1);
//            if (length < 0) {
//                break;
//            }
//            totalLength += length;
//            string name("noName");
//            string comment("noComment");
//            if (ks1->name.s != NULL) name = ks1->name.s;
//            if (ks1->comment.s != NULL) comment = ks1->comment.s;
//
//            uint64_t tuple = 0LLU, rvs_tuple = 0LLU, uni_tuple, dr_tuple, pfilter;
//            int base = 1;
//            for (int i = 0; i < length; i++) {
//                char ch = ks1->seq.s[i];
//                int basenum = BaseMap[(int)ch];
//                if (basenum != -1) {
//                    tuple = ((tuple << 2) | basenum) & tupmask;
//                    rvs_tuple = (rvs_tuple >> 2) + (((uint64_t)basenum ^ 3LLU) << rev_add_move);
//                    base++;
//                } else {
//                    base = 1;
//                }
//                if (base > kmer_size) {
//                    uni_tuple = tuple < rvs_tuple ? tuple : rvs_tuple;
//                    int dim_id = (uni_tuple & domask) >> (half_outctx_len * 2);
//
//                    if (shuffled_map.count(dim_id) == 0) {
//                        continue;
//                    }
//                    pfilter = shuffled_map[dim_id];
//                    pfilter -= dim_start;
//                    dr_tuple = (((uni_tuple & undomask0) | ((uni_tuple & undomask1) << (kmer_size * 2 - half_outctx_len * 4))) >> (drlevel * 4)) | pfilter;
//
//                    if (use64)
//                        hashValueSet64.insert(dr_tuple);
//                    else
//                        hashValueSet.insert(dr_tuple);
//                }
//            }
//        }
//
//        vector<uint32_t> hashArr;
//        vector<uint64_t> hashArr64;
//        if (use64) {
//            for (auto x : hashValueSet64) {
//                hashArr64.push_back(x);
//            }
//            tmpSketch.hashSet64 = hashArr64;
//        } else {
//            for (auto x : hashValueSet) {
//                hashArr.push_back(x);
//            }
//            tmpSketch.hashSet = hashArr;
//        }
//
//        tmpSketch.id = t;
//        gzclose(fp1);
//        kseq_destroy(ks1);
//
//        #pragma omp critical
//        {
//            sketches.push_back(tmpSketch);
//            if (t % process_bar_size == 0) {
//                cerr << "finished sketching: " << t << " genomes" << endl;
//            }
//        }
//    }
//
//    //std::sort(sketches.begin(), sketches.end(), cmpSketch);
//
//    if (!isSketchFile(outputFile)) {
//        outputFile = outputFile + ".sketch";
//    }
//
//    info.half_k = half_k;
//    info.half_subk = half_subk;
//    info.drlevel = drlevel;
//    info.id = (half_k << 8) + (half_subk << 4) + drlevel;
//    info.genomeNumber = sketches.size();
//    saveSketches(sketches, info, outputFile);
//    cerr << "save the sketches into: " << outputFile << endl;
//
//    if (!isQuery) {
//        double tstart = get_sec();
//        string dictFile = outputFile + ".dict";
//        string indexFile = outputFile + ".index";
//        transSketches(sketches, info, dictFile, indexFile, numThreads);
//        double tend = get_sec();
//        cerr << "===============the time of transSketches is: " << tend - tstart << endl;
//    }
//
//    return true;
//}
//
	
	const int Kssd::BaseMap[128] =
{
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, 0, -1, 1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, 0, -1, 1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
};
	
	kssd_parameter_t Kssd::initParameter(int half_k, int half_subk, int drlevel, int * shuffled_dim){
		// init kssd parameters
		if(half_subk - drlevel < 3){
			err(errno, "the half_subk -drlevel should at least 3\n current half_subk and drlevel aare: %d and %d respectively", half_subk, drlevel);
		}
		kssd_parameter_t parameter;
		parameter.half_k = half_k;
		parameter.half_subk = half_subk;
		parameter.drlevel = drlevel;
		int half_outctx_len = half_k - half_subk;
		parameter.half_outctx_len = half_outctx_len;
		parameter.rev_add_move = 4 * half_k - 2;
		parameter.kmer_size = 2 * half_k;

		int dim_end = 1 << 4 * (half_subk - drlevel);
		parameter.dim_start = 0;
		//parameter.dim_end = MIN_SUBCTX_DIM_SMP_SZ;
		parameter.dim_end = dim_end;
		//unsigned int hashSize = get_hashSize(half_k, drlevel);
		unsigned int hashSize = 2000;
		parameter.hashSize = hashSize;
		parameter.hashLimit = hashSize * LD_FCTR;
		parameter.shuffled_dim = shuffled_dim;

		//int component_num = half_k - drlevel > COMPONENT_SZ ? 1LU << 4 * (half_k - drlevel - COMPONENT_SZ) : 1;
		int comp_bittl = 64 - 4 * half_k;
		//cout << "the comp_bittl is: " << comp_bittl << endl;

		uint64_t tupmask = _64MASK >> comp_bittl;
		uint64_t domask = (tupmask >> (4 * half_outctx_len)) << (2 * half_outctx_len);
		uint64_t undomask = (tupmask ^ domask) & tupmask;
		uint64_t undomask1 = undomask &	(tupmask >> ((half_k + half_subk) * 2));
		uint64_t undomask0 = undomask ^ undomask1;
		parameter.domask = domask;
		parameter.tupmask = tupmask;
		parameter.undomask1 = undomask1;
		parameter.undomask0 = undomask0;
		int shuffle_dim_size = 1 << (4 * half_subk);
		//printf("the tupmask is: %lx\n", tupmask);
		//printf("the domask is: %lx\n", domask);
		//printf("the undomask0 is: %lx\n", undomask0);
		//printf("the undomask1 is: %lx\n", undomask1);

		return parameter;
	}

Kssd::Kssd(int half_k, int half_subk, int drlevel, std::unique_ptr<int[]>&& shuffled_dim)
    : half_k_(half_k), half_subk_(half_subk), drlevel_(drlevel), shuffled_dim_(std::move(shuffled_dim)),
      dim_size_(1 << (4 * half_subk)), use64_((half_k - drlevel) > 8) {
    // 初始化 shuffled_map
    for (int t = 0; t < dim_size_; t++) {
        if (shuffled_dim_[t] < (1 << (4 * (half_subk - drlevel))) && shuffled_dim_[t] >= 0) { // 确保 dim_start = 0
            shuffled_map[t] = shuffled_dim_[t];
        }
    }
}

//Kssd::Kssd(int half_k, int half_subk, int drlevel, int* shuffled_dim)
//    : half_k_(half_k), half_subk_(half_subk), drlevel_(drlevel), shuffled_dim_(shuffled_dim),
//      dim_size_(1 << (4 * half_subk)), use64_((half_k - drlevel) > 8) {
//    // 初始化 shuffled_map
//    for (int t = 0; t < dim_size_; t++) {
//        if (shuffled_dim[t] < (1 << (4 * (half_subk - drlevel))) && shuffled_dim[t] >= 0) { // 确保 dim_start = 0
//            shuffled_map[t] = shuffled_dim[t];
//        }
//    }
//}



	//Kssd::Kssd(int half_k, int half_subk, int drlevel, int* shuffled_dim, int dim_size) {
	//    // 初始化参数
	//    this->half_k = half_k;
	//    this->half_subk = half_subk;
	//    this->drlevel = drlevel;
	//    this->shuffled_dim = shuffled_dim;
	//    this->dim_size = dim_size;
	//
	//    // 调用 initParameter 来初始化 kssd_parameter_t
	//    parameter = initParameter(half_k, half_subk, drlevel, shuffled_dim);
	//
	//    // 填充 shuffled_map
	//    for (int t = 0; t < dim_size; t++) {
	//        if (shuffled_dim[t] < parameter.dim_end && shuffled_dim[t] >= parameter.dim_start) {
	//            shuffled_map.insert({t, shuffled_dim[t]});
	//            std::cout << "Inserted into shuffled_map: " << t << " -> " << shuffled_dim[t] << std::endl;
	//        }
	//    }
	//}

	// 封装的 update 方法
	//	void Kssd::update(char* seq) {
	//		uint64_t tuple = 0LLU, rvs_tuple = 0LLU, uni_tuple, dr_tuple, pfilter;
	//		int base = 1;
	//
	//		// 处理 ks1->seq.s 中的每个字符
	//		for (int i = 0; i < strlen(seq); i++) {
	//			char ch = seq[i];
	//			int basenum = BaseMap[(int)ch];
	//
	//			// 检查字符是否有效
	//			if (basenum != -1) {
	//				tuple = ((tuple << 2) | basenum) & parameter.tupmask;
	//				rvs_tuple = (rvs_tuple >> 2) + (((uint64_t)basenum ^ 3LLU) << parameter.rev_add_move);
	//				base++;
	//			} else {
	//				base = 1;
	//			}
	//
	//			if (base > parameter.kmer_size) {
	//				uni_tuple = tuple < rvs_tuple ? tuple : rvs_tuple;
	//				int dim_id = (uni_tuple & parameter.domask) >> (parameter.half_outctx_len * 2);
	//
	//				// 如果 dim_id 在 shuffled_map 中
	//				if (shuffled_map.count(dim_id) == 0) {
	//					continue;
	//				}
	//
	//				// 获取 pfilter 并更新 dr_tuple
	//				pfilter = shuffled_map[dim_id];
	//				pfilter -= parameter.dim_start;
	//				dr_tuple = (((uni_tuple & parameter.undomask0) |
	//							((uni_tuple & parameter.undomask1) << (parameter.kmer_size * 2 - parameter.half_outctx_len * 4))) >>
	//						(parameter.drlevel * 4)) |
	//					pfilter;
	//				// 其他处理逻辑，如插入哈希表等
	//				// std::cout << "dr_tuple: " << dr_tuple << std::endl;
	//				if (use64) {
	//					hashValueSet64.insert(dr_tuple);
	//				} else {
	//					hashValueSet.insert(dr_tuple);
	//				}
	//
	//			}
	//		}
	//	}

	void Kssd::update(
			const char* seq,
			int length,
			const kssd_parameter_t& param,
			bool use64,
			int kmer_size,
			std::unordered_set<uint32_t>& hashValueSet,
			std::unordered_set<uint64_t>& hashValueSet64
			) {
		uint64_t tuple = 0LLU, rvs_tuple = 0LLU, uni_tuple, dr_tuple, pfilter;
		int base = 1;

		for (int i = 0; i < length; i++) {
			char ch = seq[i];
			int basenum = (ch < 128) ? BaseMap[(int)ch] : -1; // 防止数组越界
			if (basenum != -1) {
				tuple = ((tuple << 2) | basenum) & param.tupmask;
				rvs_tuple = (rvs_tuple >> 2) + (((uint64_t)basenum ^ 3LLU) << param.rev_add_move);
				base++;
			} else {
				base = 1;
			}

			if (base > kmer_size) {
				uni_tuple = (tuple < rvs_tuple) ? tuple : rvs_tuple;
				int dim_id = (uni_tuple & param.domask) >> (param.half_outctx_len * 2);

				// 检查 dim_id 是否在 shuffled_map 中
				auto it = shuffled_map.find(dim_id);
				if (it == shuffled_map.end()) {
					continue;
				}
				pfilter = it->second;
				pfilter -= param.dim_start;

				dr_tuple = (((uni_tuple & param.undomask0) |
							((uni_tuple & param.undomask1) << (kmer_size * 2 - param.half_outctx_len * 4))) >>
						(param.drlevel * 4)) | pfilter;

				if (use64)
					hashValueSet64.insert(dr_tuple);
				else
					hashValueSet.insert(static_cast<uint32_t>(dr_tuple));
			}
		}
	}

}
void index_tridist(vector<Sketch::sketch_t>& sketches, Sketch::sketchInfo_t& info, string refSketchOut, string outputFile, int kmer_size, double maxDist, int isContainment, int numThreads){

#ifdef Timer
	double t0 = get_sec();
#endif
	string indexFile = refSketchOut + ".index";
	string dictFile = refSketchOut + ".dict";
	bool use64 = info.half_k - info.drlevel > 8 ? true : false;
	robin_hood::unordered_map<uint64_t, vector<uint32_t>> hash_map_arr;
	uint32_t* sketchSizeArr = NULL;
	size_t* offset = NULL;
	uint32_t* indexArr = NULL;
	//uint64_t* dict;
	if(use64)
	{
		cerr << "-----use hash64 in index_tridist() " << endl;
		size_t hash_number;
		FILE* fp_index = fopen(indexFile.c_str(), "rb");
		if(!fp_index){
			cerr << "ERROR: index_tridist(), cannot open index file: " << indexFile << endl;
			exit(1);
		}
		int read_hash_num = fread(&hash_number, sizeof(size_t), 1, fp_index);
		uint64_t * hash_arr = new uint64_t[hash_number];
		uint32_t * hash_size_arr = new uint32_t[hash_number];
		size_t read_hash_arr = fread(hash_arr, sizeof(uint64_t), hash_number, fp_index);
		size_t read_hash_size_arr = fread(hash_size_arr, sizeof(uint32_t), hash_number, fp_index);
		if(read_hash_num != 1 || read_hash_arr != hash_number || read_hash_size_arr != hash_number){
			cerr << "ERROR: index_tridist(), error read hash_number, hash_arr, and hash_size_arr" << endl;
			exit(1);
		}

		fclose(fp_index);

		FILE* fp_dict = fopen(dictFile.c_str(), "rb");
		if(!fp_dict){
			cerr << "ERROR: index_tridist(), cannot open dict file: " << dictFile << endl;
			exit(1);
		}
		uint32_t max_hash_size = 1LLU << 24;
		uint32_t* cur_point = new uint32_t[max_hash_size];
		for(size_t i = 0; i < hash_number; i++){
			uint32_t cur_hash_size = hash_size_arr[i];
			if(cur_hash_size > max_hash_size){
				max_hash_size = cur_hash_size;
				cur_point = new uint32_t[max_hash_size];
			}
			uint32_t hash_size = fread(cur_point, sizeof(uint32_t), cur_hash_size, fp_dict);
			if(hash_size != cur_hash_size){
				cerr << "ERROR: index_tridist(), the read hash number is not equal to the saved hash number information" << endl;
				exit(1);
			}
			vector<uint32_t> cur_genome_arr(cur_point, cur_point + cur_hash_size);
			uint64_t cur_hash = hash_arr[i];
			hash_map_arr.insert({cur_hash, cur_genome_arr});
		}
		delete [] cur_point;
		delete [] hash_arr;
		delete [] hash_size_arr;
		fclose(fp_dict);
	}
	else
	{
		cerr << "-----not use hash64 in index_tridist() " << endl;
		size_t hashSize;
		uint64_t totalIndex;
		FILE * fp_index = fopen(indexFile.c_str(), "rb");
		if(!fp_index){
			cerr << "ERROR: index_tridist(), cannot open the index sketch file: " << indexFile << endl;
			exit(1);
		}
		int read_hash_size = fread(&hashSize, sizeof(size_t), 1, fp_index);
		int read_total_index = fread(&totalIndex, sizeof(uint64_t), 1, fp_index);
		//sketchSizeArr = (uint32_t*)malloc(hashSize * sizeof(uint32_t));
		sketchSizeArr = new uint32_t[hashSize];
		size_t read_sketch_size_arr = fread(sketchSizeArr, sizeof(uint32_t), hashSize, fp_index);

		//offset = (size_t*)malloc(hashSize * sizeof(size_t));
		offset = new size_t[hashSize];
		uint64_t totalHashNumber = 0;
		for(size_t i = 0; i < hashSize; i++){
			totalHashNumber += sketchSizeArr[i];
			offset[i] = sketchSizeArr[i];
			if(i > 0) offset[i] += offset[i-1];
		}
		if(totalHashNumber != totalIndex){
			cerr << "ERROR: index_tridist(), mismatched total hash number" << endl;
			exit(1);
		}
		fclose(fp_index);

		//cerr << "the hashSize is: " << hashSize << endl;
		//cerr << "totalIndex is: " << totalIndex << endl;
		//cerr << "totalHashNumber is: " << totalHashNumber << endl;
		//cerr << "offset[n-1] is: " << offset[hashSize-1] << endl;;

		//indexArr = (uint32_t*)malloc(totalHashNumber * sizeof(uint32_t));
		indexArr = new uint32_t[totalHashNumber];
		FILE * fp_dict = fopen(dictFile.c_str(), "rb");
		if(!fp_dict){
			cerr << "ERROR: index_tridist(), cannot open the dictionary sketch file: " << dictFile << endl;
			exit(1);
		}
		size_t read_index_arr = fread(indexArr, sizeof(uint32_t), totalHashNumber, fp_dict);
		if(read_hash_size != 1 || read_total_index != 1 || read_sketch_size_arr != hashSize || read_index_arr != totalHashNumber){
			cerr << "ERROR: index_tridist(), error read hash_size, total_index, sketch_size_arr, index_arr" << endl;
			exit(1);
		}
	}

#ifdef Timer
	double t1 = get_sec();
	cerr << "===================time of read index and offset sketch file is: " << t1 - t0 << endl;
#endif
	size_t numRef = sketches.size();

	vector<FILE*> fpArr;
	vector<FILE*> fpIndexArr;
	vector<string> dist_file_list;
	vector<string> dist_index_list;
	//vector<int*> intersectionArr;
	int** intersectionArr = new int*[numThreads];

	string folderPath = outputFile + ".dir";
	string command0 = "mkdir -p " + folderPath;
	int status = system(command0.c_str());
	if(!status){
		cerr << "success create: " << folderPath << endl;
	}

	for(int i = 0; i < numThreads; i++)
	{
		string tmpName = folderPath + '/' + outputFile + '.' + to_string(i);
		dist_file_list.push_back(tmpName);
		FILE * fp0 = fopen(tmpName.c_str(), "w+");
		fpArr.push_back(fp0);

		string tmpIndexName = outputFile + ".index." + to_string(i);
		dist_index_list.push_back(tmpIndexName);
		FILE * fp1 = fopen(tmpIndexName.c_str(), "w+");
		fpIndexArr.push_back(fp1);

		//int * arr = (int*)malloc(numRef * sizeof(int));
		//int* arr = new int[numRef];
		//intersectionArr.push_back(arr);
		intersectionArr[i] = new int[numRef];
	}

	//cerr << "before generate the intersection " << endl;

	int progress_bar_size = get_progress_bar_size(numRef);
	cerr << "=====total: " << numRef << endl;
#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
	for(size_t i = 0; i < numRef; i++){
		if(i % progress_bar_size == 0) cerr << "=====finish: " << i << endl;
		int tid = omp_get_thread_num();
		fprintf(fpIndexArr[tid], "%s\t%s\n", sketches[i].fileName.c_str(), dist_file_list[tid].c_str());
		memset(intersectionArr[tid], 0, numRef * sizeof(int));
		if(use64){
			for(size_t j = 0; j < sketches[i].hashSet64.size(); j++){
				uint64_t hash64 = sketches[i].hashSet64[j];
				//if(!(dict[hash64/64] & (0x8000000000000000LLU >> (hash64 % 64))))	continue;
				if(hash_map_arr.count(hash64) == 0) continue;
				//for(auto x : hash_map_arr[hash64])
				for(size_t k = 0; k < hash_map_arr[hash64].size(); k++){
					size_t cur_index = hash_map_arr[hash64][k];
					intersectionArr[tid][cur_index]++;
					//cerr << hash64 << '\t' << cur_index << endl;
				}
			}
		}
		else{
			for(size_t j = 0; j < sketches[i].hashSet.size(); j++){
				uint32_t hash = sketches[i].hashSet[j];
				if(sketchSizeArr[hash] == 0) continue;
				size_t start = hash > 0 ? offset[hash-1] : 0;
				size_t end = offset[hash];
				for(size_t k = start; k < end; k++){
					size_t curIndex = indexArr[k];
					intersectionArr[tid][curIndex]++;
				}
			}
		}

		string strBuf("");
		for(size_t j = i+1; j < numRef; j++){
			int common = intersectionArr[tid][j];
			int size0, size1;
			if(use64){
				size0 = sketches[i].hashSet64.size();
				size1 = sketches[j].hashSet64.size();
			}
			else{
				size0 = sketches[i].hashSet.size();
				size1 = sketches[j].hashSet.size();
			}
			if(!isContainment){
				int denom = size0 + size1 - common;
				double jaccard;
				if(size0 == 0 || size1 == 0)
					jaccard = 0.0;
				else
					jaccard = (double)common / denom;
				double mashD;
				if(jaccard == 1.0)
					mashD = 0.0;
				else if(jaccard == 0.0)
					mashD = 1.0;
				else
					mashD = (double)-1.0 / kmer_size * log((2 * jaccard)/(1.0 + jaccard));
				if(mashD < maxDist){
					strBuf += sketches[j].fileName + '\t' + sketches[i].fileName + '\t' + to_string(common) + '|' + to_string(size0) + '|' + to_string(size1) + '\t' + to_string(jaccard) + '\t' + to_string(mashD) + '\n';
				}
				//fprintf(fpArr[tid], " %s\t%s\t%d|%d|%d\t%lf\t%lf\n", sketches[j].fileName.c_str(), sketches[i].fileName.c_str(), common, size0, size1, jaccard, mashD);
			}
			else{
				int denom = std::min(size0, size1);
				double containment;
				if(size0 == 0 || size1 == 0)
					containment = 0.0;
				else
					containment = (double)common / denom;
				double AafD;
				if(containment == 1.0)
					AafD = 0.0;
				else if(containment == 0.0)
					AafD = 1.0;
				else
					AafD = (double)-1.0 / kmer_size * log(containment);
				if(AafD < maxDist)
					strBuf += sketches[j].fileName + '\t' + sketches[i].fileName + '\t' + to_string(common) + '|' + to_string(size0) + '|' + to_string(size1) + '\t' + to_string(containment) + '\t' + to_string(AafD) + '\n';
				//fprintf(fpArr[tid], " %s\t%s\t%d|%d|%d\t%lf\t%lf\n", sketches[j].fileName.c_str(), sketches[i].fileName.c_str(), common, size0, size1, containment, AafD);
			}
		}
		fprintf(fpArr[tid], "%s", strBuf.c_str());
		strBuf = "";
	}


	//cerr << "finished multithread computing" << endl;

	for(int i = 0; i < numThreads; i++)
	{
		fclose(fpArr[i]);
		fclose(fpIndexArr[i]);
		delete [] intersectionArr[i];
	}
	//cerr << "finished fpArr fclose" << endl;

#ifdef Timer
	double t2 = get_sec();
	cerr << "===================time of multiple threads distance computing and save the subFile is: " << t2 - t1 << endl;
#endif

	uint64_t totalSize = 0;
	uint64_t maxSize = 1LLU << 32; //max total distance file size 4GB
	bool isMerge = false;
	for(int i = 0; i < numThreads; i++){
		struct stat cur_stat;
		stat(dist_file_list[i].c_str(), &cur_stat);
		uint64_t curSize = cur_stat.st_size;
		totalSize += curSize;
	}
	if(totalSize <= maxSize)	isMerge = true;

	if(isMerge){
		FILE * cofp;
		FILE * com_cofp = fopen(outputFile.c_str(), "w");
		cerr << "-----save the output distance file: " << outputFile << endl;
		fprintf(com_cofp, " genome0\tgenome1\tcommon|size0|size1\tjaccard\tmashD\n");
		int bufSize = 1 << 24;
		int lengthRead = 0;
		char * bufRead = (char*)malloc((bufSize+1) * sizeof(char));
		for(int i = 0; i < numThreads; i++)
		{
			cofp = fopen(dist_file_list[i].c_str(), "rb+");
			while(1)
			{
				lengthRead = fread(bufRead, sizeof(char), bufSize, cofp);
				//cerr << "the lengthRead is: " << lengthRead << endl;
				fwrite(bufRead, sizeof(char), lengthRead, com_cofp);
				if(lengthRead < bufSize) break;
			}
			fclose(cofp);
			remove(dist_file_list[i].c_str());
			remove(dist_index_list[i].c_str());
		}
		remove(folderPath.c_str());

		free(bufRead);
		fclose(com_cofp);
	}
	else{
		cerr << "-----the output distance file is too big to merge into one single file, saving the result into directory: " << folderPath << endl;
		FILE * cofp1;
		string outputIndexFile = outputFile + ".index";
		cerr << "-----save the index between genomes and distance sub-files into: " << outputIndexFile << endl;
		FILE * com_cofp1 = fopen(outputIndexFile.c_str(), "w+");
		fprintf(com_cofp1, "genomeName\tdistFileName\n");
		int bufSize = 1 << 24;
		int lengthRead = 0;
		char * bufRead = (char*)malloc((bufSize+1) * sizeof(char));
		for(int i = 0; i < numThreads; i++){
			cofp1 = fopen(dist_index_list[i].c_str(), "rb+");
			while(1){
				lengthRead = fread(bufRead, sizeof(char), bufSize, cofp1);
				fwrite(bufRead, sizeof(char), lengthRead, com_cofp1);
				if(lengthRead < bufSize) break;
			}
			fclose(cofp1);
			remove(dist_index_list[i].c_str());
		}
		free(bufRead);
		fclose(com_cofp1);
	}

#ifdef Timer
	double t3 = get_sec();
	cerr << "===================time of merge the subFiles into final files is: " << t3 - t2 << endl;
#endif

}

