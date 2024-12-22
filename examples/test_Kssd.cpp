#include <iostream>
#include <errno.h>
#include <err.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <string>
#include <vector>
#include <stdint.h>
#include <unordered_set>
#include <unordered_map>
#include "kseq.h"
#include "zlib.h"
#include "shuffle.h"
#include "Kssd.cpp"
#include "Sketch.h"
#include "robin_hood.h"
#include <algorithm>
#include <omp.h>
#include <cstdlib>
#include <sys/stat.h> 
KSEQ_INIT(gzFile, gzread);

using namespace Sketch;

bool cmpFile(fileInfo_t f1, fileInfo_t f2){
	return f1.fileSize > f2.fileSize;
}

int main() {
	//default
	int half_k = 10;
	int half_subk = 6;
	int drlevel = 3;
	sketchInfo_t info;
	auto result = Sketch::read_shuffled_file("/home/user_home/zt/bioRabbitSketch/RabbitSketch/src/shuf_file/L3K10.shuf");
	half_k = std::get<0>(result);
	half_subk = std::get<1>(result);
	drlevel = std::get<2>(result);
	int* shuffled_dim_ptr = std::get<3>(result);
	//std::unique_ptr<int[]> shuffled_dim_ptr = std::move(std::get<3>(result));
	if (half_k == -1) {
		std::cerr << "Error reading shuffled file" << std::endl;
		return -1;
	}
	std::cout << "half_k: " << half_k << std::endl;
	std::cout << "half_subk: " << half_subk << std::endl;
	std::cout << "drlevel: " << drlevel << std::endl;


	kssd_parameter_t kssdPara(half_k, half_subk, drlevel, shuffled_dim_ptr);
	vector<Sketch::Kssd *> vkssd;	
	bool isQuery = false;
	std::ifstream fs("100.list");
	if (!fs) {
		std::cerr << "cant open 100.list" << std::endl;
		return -1;
	}

	std::vector<fileInfo_t> fileList;
	uint64_t totalSize = 0;
	std::string fileName;
	while (std::getline(fs, fileName)) {
		struct stat cur_stat;
		if (stat(fileName.c_str(), &cur_stat) != 0) {
			std::cerr << "err" << fileName << std::endl;
			continue;
		}
		uint64_t curSize = cur_stat.st_size;
		totalSize += curSize;
		fileInfo_t tmpF;
		tmpF.fileName = fileName;
		tmpF.fileSize = curSize;
		fileList.push_back(tmpF);
	}
	fs.close();
	std::sort(fileList.begin(), fileList.end(), cmpFile);
	int small_file_number = fileList.size();
	int process_bar_size = get_progress_bar_size(small_file_number); 
	std::cerr << "===== total file num :" << small_file_number << std::endl;
#pragma omp parallel for schedule(dynamic)
	for (int t = 0; t < small_file_number; t++) {
		Sketch::Kssd* kssd = new Sketch::Kssd(kssdPara);
		Sketch::sketch_t tmpSketch;
		kssd->fileName = fileList[t].fileName;
		kssd->id = t;
		//tmpSketch.fileName = fileList[t].fileName;
		//tmpSketch.id = t;

		gzFile fp1 = gzopen(fileList[t].fileName.c_str(), "r");
		kseq_t* ks1 = kseq_init(fp1);
		while (1) {
			int length = kseq_read(ks1);
			if (length < 0) {
				break;
			}

			kssd->update(ks1->seq.s);
		}
#pragma omp critical
		{
			vkssd.push_back(kssd);
			if (t % process_bar_size == 0) {
				std::cerr << "finish sketch: " << t << " genome" << std::endl;
			}
		}


		gzclose(fp1);
		kseq_destroy(ks1);
	}

	if (!isSketchFile("100.sketch")) { 
		std::cerr << "err:cant creat file" << std::endl;
		return -1;
	}

	//	for(int i=0;i<100;i++){
	//		std::cout << vkssd[i]->distance(vkssd[i])<< "aaa"<<endl;
	//	}
	info.half_k = half_k;
	info.half_subk = half_subk;
	info.drlevel = drlevel;
	std::string outputFile = "100.sketch";
	info.id = (half_k << 8) + (half_subk << 4) + drlevel;
	//	info.genomeNumber = sketches.size();
	std::cout << "sketch num: " << vkssd.size() << std::endl;

	saveSketches(vkssd, info, "100.sketch"); 
	std::cerr << "save sketches to : 100.sketch" << std::endl;

	if (!isQuery) {
		double tstart = get_sec();
		std::string dictFile = outputFile + ".dict";
		std::string indexFile = outputFile + ".index";
		transSketches(vkssd, info, dictFile, indexFile, 64); 
		double tend = get_sec();
		std::cerr << "=============== transSketches time: " << tend - tstart << " s" << std::endl;
	}
	//			delete[] shuffled_info->shuffled_dim;
	//			delete shuffled_info;
	Sketch::index_tridist(vkssd, info, "100.sketch", "100.sketch.dist", 20, 1.0, 0, 64);
	return 0;
}


