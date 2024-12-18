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
        int half_k = 10;
        int half_subk = 6;
        int drlevel = 3;
        std::vector<Sketch::sketch_t> sketches;
        sketchInfo_t info;
        dim_shuffle_t* shuffled_info = read_shuffle_dim("/home/user_home/zt/bioRabbitSketch/RabbitSketch/src/shuf_file/L3K10.shuf");
        if (!shuffled_info) {
            std::cerr << "无法读取shuffle文件。" << std::endl;
            return -1;
        }
        half_k = shuffled_info->dim_shuffle_stat.k;
        half_subk = shuffled_info->dim_shuffle_stat.subk;
        drlevel = shuffled_info->dim_shuffle_stat.drlevel;
        std::unique_ptr<int[]> shuffled_dim_ptr(new int[1 << (4 * half_subk)]);
        std::copy(shuffled_info->shuffled_dim, shuffled_info->shuffled_dim + (1 << (4 * half_subk)), shuffled_dim_ptr.get());
        Sketch::Kssd kssd(half_k, half_subk, drlevel, std::move(shuffled_dim_ptr));
				Sketch::kssd_parameter_t kssd_parameter = kssd.initParameter(half_k, half_subk, drlevel, shuffled_dim_ptr.get());
        int half_outctx_len = kssd_parameter.half_outctx_len;
        int* shuffled_dim = kssd_parameter.shuffled_dim;
        int kmer_size = kssd_parameter.kmer_size;
        bool use64 = (half_k - drlevel) > 8;
        bool isQuery = false;
        int dim_size = 1 << (4 * (half_k - half_outctx_len));
        std::ifstream fs("100.list");
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
        std::sort(fileList.begin(), fileList.end(), cmpFile);
        int small_file_number = fileList.size();
        int process_bar_size = get_progress_bar_size(small_file_number); 
        std::cerr << "===== 总小文件数量: " << small_file_number << std::endl;
        #pragma omp parallel for schedule(dynamic)
        for (int t = 0; t < small_file_number; t++) {
            Sketch::sketch_t tmpSketch;
            tmpSketch.fileName = fileList[t].fileName;
            tmpSketch.id = t;

            gzFile fp1 = gzopen(fileList[t].fileName.c_str(), "r");
            if (fp1 == NULL) {
                #pragma omp critical
                {
                    std::cerr << "无法打开基因组文件: " << fileList[t].fileName << std::endl;
                }
                continue; 
            }

            kseq_t* ks1 = kseq_init(fp1);

            std::unordered_set<uint32_t> hashValueSet;
            std::unordered_set<uint64_t> hashValueSet64;

            while (1) {
                int length = kseq_read(ks1);
                if (length < 0) {
                    break;
                }

                const char* seq = ks1->seq.s;
                kssd.update(seq, length, kssd_parameter, use64, kmer_size, hashValueSet, hashValueSet64);
            }

            if (use64) {
                tmpSketch.hashSet64.assign(hashValueSet64.begin(), hashValueSet64.end());
            } else {
                tmpSketch.hashSet.assign(hashValueSet.begin(), hashValueSet.end());
            }

            gzclose(fp1);
            kseq_destroy(ks1);

            #pragma omp critical
            {
                sketches.push_back(tmpSketch);
                if (t % process_bar_size == 0) {
                    std::cerr << "已完成sketch生成: " << t << " 个基因组" << std::endl;
                }
            }
        }

        if (!isSketchFile("100.sketch")) { 
            std::cerr << "错误: 无法创建 sketch 文件。" << std::endl;
            return -1;
        }

        info.half_k = half_k;
        info.half_subk = half_subk;
        info.drlevel = drlevel;
        std::string outputFile = "100.sketch";
        info.id = (half_k << 8) + (half_subk << 4) + drlevel;
        info.genomeNumber = sketches.size();
        std::cout << "sketch 数量: " << sketches.size() << std::endl;

        saveSketches(sketches, info, "100.sketch"); 
        std::cerr << "已保存 sketches 到: 100.sketch" << std::endl;

        if (!isQuery) {
            double tstart = get_sec();
            std::string dictFile = outputFile + ".dict";
            std::string indexFile = outputFile + ".index";
            transSketches(sketches, info, dictFile, indexFile, 64); 
            double tend = get_sec();
            std::cerr << "=============== transSketches 运行时间: " << tend - tstart << " 秒" << std::endl;
        }
        delete[] shuffled_info->shuffled_dim;
        delete shuffled_info;
				index_tridist(sketches, info, "100.sketch", "100.sketch.dist", 20, 1.0, 0, 64);
        return 0;
    }


