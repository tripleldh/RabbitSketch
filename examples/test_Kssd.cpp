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
        // 声明变量
        int half_k = 10;
        int half_subk = 6;
        int drlevel = 3;
        std::vector<Sketch::sketch_t> sketches;
        sketchInfo_t info;

        // 调用读取 shuffle 文件的函数
        dim_shuffle_t* shuffled_info = read_shuffle_dim("/home/user_home/zt/bioRabbitSketch/RabbitSketch/src/shuf_file/L3K10.shuf");
        if (!shuffled_info) {
            std::cerr << "无法读取shuffle文件。" << std::endl;
            return -1;
        }

        // 使用读取到的值来设置参数
        half_k = shuffled_info->dim_shuffle_stat.k;
        half_subk = shuffled_info->dim_shuffle_stat.subk;
        drlevel = shuffled_info->dim_shuffle_stat.drlevel;

        // 将 shuffled_dim 转换为 unique_ptr
        std::unique_ptr<int[]> shuffled_dim_ptr(new int[1 << (4 * half_subk)]);
        std::copy(shuffled_info->shuffled_dim, shuffled_info->shuffled_dim + (1 << (4 * half_subk)), shuffled_dim_ptr.get());

        // 初始化 Kssd 参数
        Sketch::Kssd kssd(half_k, half_subk, drlevel, std::move(shuffled_dim_ptr));
				Sketch::kssd_parameter_t kssd_parameter = kssd.initParameter(half_k, half_subk, drlevel, shuffled_dim_ptr.get());

        // 生成 sketch 文件参数
        int rev_add_move = kssd_parameter.rev_add_move;
        int half_outctx_len = kssd_parameter.half_outctx_len;
        int* shuffled_dim = kssd_parameter.shuffled_dim;
        int dim_start = kssd_parameter.dim_start;
        int dim_end = kssd_parameter.dim_end;
        int kmer_size = kssd_parameter.kmer_size;
        uint64_t tupmask = kssd_parameter.tupmask;
        uint64_t undomask0 = kssd_parameter.undomask0;
        uint64_t undomask1 = kssd_parameter.undomask1;
        uint64_t domask = kssd_parameter.domask;

        bool use64 = (half_k - drlevel) > 8;
        bool isQuery = false;
        int dim_size = 1 << (4 * (half_k - half_outctx_len));

        // 打开文件列表
        std::ifstream fs("bacteria.list");
        if (!fs) {
            std::cerr << "无法打开输入文件: 100.list" << std::endl;
            return -1;
        }

        std::vector<fileInfo_t> fileList;
        uint64_t totalSize = 0;
        std::string fileName;
        while (std::getline(fs, fileName)) {
            struct stat cur_stat;
            if (stat(fileName.c_str(), &cur_stat) != 0) {
                std::cerr << "无法获取文件状态: " << fileName << std::endl;
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

        // 按文件大小排序
        std::sort(fileList.begin(), fileList.end(), cmpFile);

        // 所有文件现在都被视为小文件，无需区分大/小文件
        int small_file_number = fileList.size();
        int process_bar_size = get_progress_bar_size(small_file_number); // 确保这个函数已定义
        std::cerr << "===== 总小文件数量: " << small_file_number << std::endl;

        // 使用 OpenMP 并行处理每个文件
        #pragma omp parallel for schedule(dynamic)
        for (int t = 0; t < small_file_number; t++) {
            Sketch::sketch_t tmpSketch;
            tmpSketch.fileName = fileList[t].fileName;
            tmpSketch.id = t;

            // 打开压缩的FASTA文件
            gzFile fp1 = gzopen(fileList[t].fileName.c_str(), "r");
            if (fp1 == NULL) {
                #pragma omp critical
                {
                    std::cerr << "无法打开基因组文件: " << fileList[t].fileName << std::endl;
                }
                continue; // 继续处理下一个文件
            }

            kseq_t* ks1 = kseq_init(fp1);

            // 创建每个线程独立的哈希集合
            std::unordered_set<uint32_t> hashValueSet;
            std::unordered_set<uint64_t> hashValueSet64;

            while (1) {
                int length = kseq_read(ks1);
                if (length < 0) {
                    break;
                }

                const char* seq = ks1->seq.s;
                // 调用 Kssd::update 函数处理序列
                kssd.update(
                    seq,
                    length,
                    kssd_parameter,
                    use64,
                    kmer_size,
                    hashValueSet,
                    hashValueSet64
                );
            }

            // 将哈希值保存到 tmpSketch 中
            if (use64) {
                tmpSketch.hashSet64.assign(hashValueSet64.begin(), hashValueSet64.end());
            } else {
                tmpSketch.hashSet.assign(hashValueSet.begin(), hashValueSet.end());
            }

            gzclose(fp1);
            kseq_destroy(ks1);

            // 保护对 sketches 向量的并发访问
            #pragma omp critical
            {
                sketches.push_back(tmpSketch);
                if (t % process_bar_size == 0) {
                    std::cerr << "已完成sketch生成: " << t << " 个基因组" << std::endl;
                }
            }
        }

        // 检查 sketch 文件是否创建成功
        if (!isSketchFile("bacteria.sketch")) { // 确保 isSketchFile 已定义
            std::cerr << "错误: 无法创建 sketch 文件。" << std::endl;
            return -1;
        }

        // 设置 sketch 信息
        info.half_k = half_k;
        info.half_subk = half_subk;
        info.drlevel = drlevel;
        std::string outputFile = "bacteria.sketch";
        info.id = (half_k << 8) + (half_subk << 4) + drlevel;
        info.genomeNumber = sketches.size();
        std::cout << "sketch 数量: " << sketches.size() << std::endl;

        // 保存 sketches
        saveSketches(sketches, info, "bacteria.sketch"); // 确保 saveSketches 已定义
        std::cerr << "已保存 sketches 到: bacteria.sketch" << std::endl;

        if (!isQuery) {
            double tstart = get_sec(); // 确保 get_sec 已定义
            std::string dictFile = outputFile + ".dict";
            std::string indexFile = outputFile + ".index";
            transSketches(sketches, info, dictFile, indexFile, 64); // 确保 transSketches 已定义
            double tend = get_sec();
            std::cerr << "=============== transSketches 运行时间: " << tend - tstart << " 秒" << std::endl;
        }

        // 生成索引分布
        //index_tridist(sketches, info, "bacteria.sketch", "bacteria.sketch.dist", 20, 1.0, 0, 64); // 确保 index_tridist 已定义

        // 释放 shuffled_info 的内存（假设 read_shuffle_dim 使用 new 分配）
        delete[] shuffled_info->shuffled_dim;
        delete shuffled_info;

        return 0;
    }

//int main() {
//	// 声明变量
//	int half_k = 10;
//	int half_subk = 6;
//	int drlevel = 3;
//	vector<sketch_t> sketches;
//	sketchInfo_t info;
//
//	// 调用读取 shuffle 文件的函数
//	dim_shuffle_t* shuffled_info = read_shuffle_dim("/home/user_home/zt/bioRabbitSketch/RabbitSketch/src/shuf_file/L3K10.shuf");
//
//	// 使用读取到的值来设置参数
//	half_k = shuffled_info->dim_shuffle_stat.k;
//	half_subk = shuffled_info->dim_shuffle_stat.subk;
//	drlevel = shuffled_info->dim_shuffle_stat.drlevel;
//	//int dim_size = 16777216;
//	//int dim_size = 1 << 4 * (half_subk - drlevel);
//	Sketch::Kssd kssd(half_k, half_subk, drlevel, shuffled_info->shuffled_dim);
//	// 初始化 kssd 参数
//	kssd_parameter_t kssd_parameter = kssd.initParameter(half_k, half_subk, drlevel, shuffled_info->shuffled_dim);
//
//	// 生成 sketch 文件
//	int rev_add_move = kssd_parameter.rev_add_move;
//	int half_outctx_len = kssd_parameter.half_outctx_len;
//	int* shuffled_dim = kssd_parameter.shuffled_dim;
//	int dim_start = kssd_parameter.dim_start;
//	int dim_end = kssd_parameter.dim_end;
//	int kmer_size = kssd_parameter.kmer_size;
//	uint64_t tupmask = kssd_parameter.tupmask;
//	uint64_t undomask0 = kssd_parameter.undomask0;
//	uint64_t undomask1 = kssd_parameter.undomask1;
//	uint64_t domask = kssd_parameter.domask;
////
////std::cout << "half_outctx_len: " << half_outctx_len << std::endl;
////std::cout << "shuffled_dim: ";
////std::cout << "dim_start: " << dim_start << std::endl;
////std::cout << "dim_end: " << dim_end << std::endl;
////std::cout << "kmer_size: " << kmer_size << std::endl;
//////std::cout << "dim_size: " << dim_size << std::endl;
////std::cout << "tupmask: " << std::hex << tupmask << std::dec << std::endl;
////std::cout << "undomask0: " << std::hex << undomask0 << std::dec << std::endl;
////std::cout << "undomask1: " << std::hex << undomask1 << std::dec << std::endl;
////std::cout << "domask: " << std::hex << domask << std::dec << std::endl;
//	bool use64 = half_k - drlevel > 8 ? true : false;
//	bool isQuery = false;
//	int dim_size = 1 << 4 * (half_k - half_outctx_len);
////	robin_hood::unordered_map<uint32_t, int> shuffled_map;
////	for (int t = 0; t < dim_size; t++) {
////		if (shuffled_dim[t] < dim_end && shuffled_dim[t] >= dim_start) {
////			shuffled_map.insert({t, shuffled_dim[t]});
////		}
////	}
//
//
//	ifstream fs("100.list");
//	if (!fs) {
//		std::cerr << "cannot open the inputFile: bacteria.list" << std::endl;
//		return -1;
//	}
//	vector<fileInfo_t> fileList;
//	uint64_t totalSize = 0;
//	string fileName;
//	while (getline(fs, fileName)) {
//		struct stat cur_stat;
//		stat(fileName.c_str(), &cur_stat);
//		uint64_t curSize = cur_stat.st_size;
//		totalSize += curSize;
//		fileInfo_t tmpF;
//		tmpF.fileName = fileName;
//		tmpF.fileSize = curSize;
//		fileList.push_back(tmpF);
//	}
//	std::sort(fileList.begin(), fileList.end(), cmpFile);
//
//	// All files are considered small files now, so no need for big/small file differentiation
//	int small_file_number = fileList.size();
//	int process_bar_size = get_progress_bar_size(small_file_number);
//	std::cerr << "=====total small files: " << small_file_number << std::endl;
//
//#pragma omp parallel for num_threads(64) schedule(dynamic)
//	for (size_t t = 0; t < small_file_number; t++) {
//		sketch_t tmpSketch;
//		gzFile fp1;
//		kseq_t* ks1;
//		fp1 = gzopen(fileList[t].fileName.c_str(), "r");
//		if (fp1 == NULL) {
//			//std::cerr << "cannot open the genome file: " << fileList[t].fileName << std::endl;
//			//continue;
//			err(errno, "cannot open the genome file: %s\n", fileList[t].fileName.c_str());
//		}
//		ks1 = kseq_init(fp1);
//		uint64_t totalLength = 0;
//		tmpSketch.fileName = fileList[t].fileName;
//
//		unordered_set<uint32_t> hashValueSet;
//		unordered_set<uint64_t> hashValueSet64;
////		kssd.hashValueSet.clear();
////		kssd.hashValueSet64.clear();
//		while (1) {
//			int length = kseq_read(ks1);
//			//std::cout << "Read length: " << length << std::endl;
//			if (length < 0) {
//				break;
//			}
//			totalLength += length;
//			string name("noName");
//			string comment("noComment");
//			if (ks1->name.s != NULL) name = ks1->name.s;
//			if (ks1->comment.s != NULL) comment = ks1->comment.s;
////			kssd.update(ks1->seq.s);
//			uint64_t tuple = 0LLU, rvs_tuple = 0LLU, uni_tuple, dr_tuple, pfilter;
//			int base = 1;
//			for (int i = 0; i < length; i++) {
//				char ch = ks1->seq.s[i];
//				int basenum = kssd.BaseMap[(int)ch];
//				if (basenum != -1) {
//					tuple = ((tuple << 2) | basenum) & tupmask;
//					rvs_tuple = (rvs_tuple >> 2) + (((uint64_t)basenum ^ 3LLU) << rev_add_move);
//					base++;
//				} else {
//					base = 1;
//				}
//				if (base > kmer_size) {
//					uni_tuple = tuple < rvs_tuple ? tuple : rvs_tuple;
//					int dim_id = (uni_tuple & domask) >> (half_outctx_len * 2);
//					//std::cout << "dim_id: " << dim_id << std::endl;
//					if (kssd.shuffled_map.count(dim_id) == 0) {
//						continue;
//					}
//					pfilter = kssd.shuffled_map[dim_id];
//					pfilter -= dim_start;
//					dr_tuple = (((uni_tuple & undomask0) | ((uni_tuple & undomask1) << (kmer_size * 2 - half_outctx_len * 4))) >> (drlevel * 4)) | pfilter;
//					//std::cout << "dr_tuple: " << dr_tuple << std::endl;
//
//					if (use64)
//						hashValueSet64.insert(dr_tuple);
//					else
//						hashValueSet.insert(dr_tuple);
//				}
//			}
//		}
//
//		vector<uint32_t> hashArr;
//		vector<uint64_t> hashArr64;
//		if (use64) {
//			for (auto x : hashValueSet64) {
//				hashArr64.push_back(x);
//			}
//			tmpSketch.hashSet64 = hashArr64;
//		} else {
//			for (auto x : hashValueSet) {
//				hashArr.push_back(x);
//			}
//			tmpSketch.hashSet = hashArr;
//		}
//
//
//		tmpSketch.id = t;
//		gzclose(fp1);
//		kseq_destroy(ks1);
//
//#pragma omp critical
//		{
//			sketches.push_back(tmpSketch);
//			if (t % process_bar_size == 0) {
//				std::cerr << "finished sketching: " << t << " genomes" << std::endl;
//			}
//		}
//	}
//
//	//std::sort(sketches.begin(), sketches.end(), cmpSketch);
//
//	if (!isSketchFile("bacteria.sketch")) {
//		std::cerr << "Error: The sketch file could not be created." << std::endl;
//		return -1;
//	}
//
//	info.half_k = half_k;
//	info.half_subk = half_subk;
//	info.drlevel = drlevel;
//	string outputFile = "bacteria.sketch";
//	info.id = (half_k << 8) + (half_subk << 4) + drlevel;
//	info.genomeNumber = sketches.size();
//	std::cout << "size: " << sketches.size() << std::endl;
//	saveSketches(sketches, info, "bacteria.sketch");
//	std::cerr << "save the sketches into: bacteria.sketch" << std::endl;
//	if(!isQuery){
//		double tstart = get_sec();
//		string dictFile = outputFile + ".dict";
//		string indexFile = outputFile + ".index";
//		transSketches(sketches, info, dictFile, indexFile, 64);
//		double tend = get_sec();
//		cerr << "===============the time of transSketches is: " << tend - tstart << endl;
//	}
//
//	index_tridist(sketches, info, "bacteria.sketch", "bacteria.sketch.dist", 20, 1.0, 0, 64);
//
//	return 0;
//}
