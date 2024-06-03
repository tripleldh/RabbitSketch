#include "Sketch.h"
//#include <iostream>
#include <sys/time.h>
#include <zlib.h>
#include "kseq.h"
#include <vector>
#include <math.h>
#include <random>
//#include <sstream>
#include <fstream>
#include <err.h>
#include <sys/stat.h>
#include <omp.h>
#include <sstream>
using namespace std;

//struct HLLInfo{
//	string name;
//	Sketch::HyperLogLog hll;
//};
typedef struct fileInfo
{
  string fileName;
  //	uint64_t fileSize;
} fileInfo_t;

KSEQ_INIT(gzFile, gzread)

  double get_sec(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + (double)tv.tv_usec/1000000;
  }


int main(int argc, char* argv[])
{
	if(argc < 3){
		cerr << "run as: " << argv[0] << " bac.txt" << endl;
		return 1;
	}
  string inputFile = argv[1];
  double thres = stod(argv[2]);
  //string inputFile = "bac.txt";
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
    //	uint64_t curSize = cur_stat.st_size;
    //	totalSize += curSize;
    fileInfo_t tmpF;
    tmpF.fileName = fileName;
    //	tmpF.fileSize = curSize;
    fileList.push_back(tmpF);
  }
  //	std::sort(fileList.begin(), fileList.end(), cmpFile);
  //	uint64_t limitSize = totalSize / numThreads;
  //cerr << "the total fileNumber is: " << fileList.size() << endl;
  //cerr << "the big fasta file is: " << numBigFasta << endl;
  //cerr << "the small fileNumber is: " << smallFileArr.size() << endl;
  double t1 = get_sec();
  vector<string> fileArr;
  for(size_t i = 0; i < fileList.size(); i++){
    fileArr.push_back(fileList[i].fileName);
  }
  //fileArr.resize(10); 
  int small_file_number = fileArr.size();
  //vector<Sketch::MinHash *> vmh;
 
	//int half_k = 10;
	//int half_subk = 6;
	//int drlevel = 3;
	//Sketch::KSSDParameters kssdPara(half_k, half_subk, drlevel);
  vector<Sketch::HyperLogLog> vhlog;
	//vector<Sketch::KSSD *> vkssd;
	//vector<HLLInfo> vhlog;
 
 //Sketch::WMHParameters parameter;
 //parameter.kmerSize = 21;
 //parameter.sketchSize = 10;
 //parameter.windowSize = 20;
 //parameter.r = (double *)malloc(parameter.sketchSize * pow(parameter.kmerSize, 4) * sizeof(double));
 //parameter.c = (double *)malloc(parameter.sketchSize * pow(parameter.kmerSize, 4) * sizeof(double));
 //parameter.b = (double *)malloc(parameter.sketchSize * pow(parameter.kmerSize, 4) * sizeof(double));
 //getCWS(parameter.r, parameter.c, parameter.b, parameter.sketchSize, pow(parameter.kmerSize, 4));
 //int half_k = 10;
 //int half_subk = 6;
 //int drlevel = 3;
 //vector<Sketch::WMinHash *>vwmh;
  //int process_bar_size = get_progress_bar_size(small_file_number);
  cerr << "=====total small files: " << small_file_number << endl;
	vector<string> resFileName;
#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
  for(size_t t = 0; t <small_file_number; t++)
  {
    //int tid = omp_get_thread_num();
    gzFile fp1;
    kseq_t * ks1;
    fp1 = gzopen(fileArr[t].c_str(), "r");
    if(fp1 == NULL){
      err(errno, "cannot open the genome file: %s\n", fileArr[t].c_str());
    }
    ks1 = kseq_init(fp1);
    
		//Sketch::MinHash * mh1 = new Sketch::MinHash();
		
    static const size_t BITS = 10; //24
    Sketch::HyperLogLog hll(BITS);
		while(1){
      int length = kseq_read(ks1);
      if(length < 0){
        break;
 			}

		//mh1->update(ks1->seq.s);	
		//Sketch::KSSD * kssd = new Sketch::KSSD(kssdPara);
		//kssd->update(ks1->seq.s);
      //Sketch::WMinHash * wmh1 = new Sketch::WMinHash(parameter);
      //wmh1->update(ks1->seq.s);
      //static const size_t BITS = 10; //24
      //Sketch::HyperLogLog hll(BITS);
     hll.update(ks1->seq.s);
#pragma omp critical
      {
        //vwmh.push_back(wmh1);
        vhlog.push_back(hll);

				//vmh.push_back(mh1);
		  //	vkssd.push_back(kssd);
				resFileName.push_back(fileArr[t]);
				//if(vhlog.size() >= 100) exit(1);
      }
      //Sketch::MinHash * mh1 = new Sketch::MinHash();	
      //mh1->update(ks1->seq.s);

      //vmh.push_back(mh1);

    }//end while, read the file
    gzclose(fp1);
    kseq_destroy(ks1);
  }
double t2 = get_sec();
cerr << "sketch time is: " << t2 - t1 << endl;

vector<FILE*> fp_arr(numThreads); // 预先分配向量
string prefixName = "zt_res/res.dist.";

#define MAX_BUFFER_SIZE 400*1024*1024

#pragma omp parallel num_threads(numThreads)
{
    int tid = omp_get_thread_num();
    string fileName = prefixName + to_string(tid);
    fp_arr[tid] = fopen(fileName.c_str(), "w+");
    
    stringstream buffer;
    
    #pragma omp for schedule(dynamic)
    for (int i = 0; i < fileArr.size(); i++) {
        for (int j = i + 1; j < fileArr.size(); j++) {
            //double distance1 = vkssd[i]->distance(vkssd[j]);
						double distance1 = vhlog[i].distance(vhlog[j]);
	    if(distance1 <thres){
            buffer << resFileName[i] << "\t" << resFileName[j] << "\t" << distance1 << "\n";
            if (buffer.tellp() >= MAX_BUFFER_SIZE || (i == fileArr.size() - 1 && j == fileArr.size() - 1)) {
                fprintf(fp_arr[tid], "%s", buffer.str().c_str());
                fflush(fp_arr[tid]);
		buffer.str(""); // 清空缓冲区
            }
        }
      }
    }
     
   if (buffer.tellp() > 0) {
	 fprintf(fp_arr[tid], "%s", buffer.str().c_str());
         fflush(fp_arr[tid]);
    }
   fclose(fp_arr[tid]);
}

double t3 = get_sec();
cerr << "dist time is: " << t3 - t2 << endl;

}




















