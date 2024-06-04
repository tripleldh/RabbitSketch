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


typedef struct fileInfo
{
  string fileName;
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
 
	int half_k = 10;
	int half_subk = 6;
	int drlevel = 3;
	Sketch::KSSDParameters kssdPara(half_k, half_subk, drlevel);
	vector<Sketch::KSSD *> vkssd;
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
    Sketch::KSSD * kssd = new Sketch::KSSD(kssdPara);
		while(1){
      int length = kseq_read(ks1);
      if(length < 0){
        break;
 			}

		kssd->update(ks1->seq.s);
		#pragma omp critical
      {
		  	vkssd.push_back(kssd);
				resFileName.push_back(fileArr[t]);
      }

    }//end while, read the file
    gzclose(fp1);
    kseq_destroy(ks1);
  }
  double t2 = get_sec();
  cerr << "sketch time is: " << t2 - t1 << endl;
  
  vector<FILE*> fp_arr(numThreads);

	string cmd = "mkdir -p res_dir";
	int ret = system(cmd.c_str());
	if(ret == 0){
		cerr << "write the result int the directory: res_dir" << endl;
	}
	else{
		cerr << "ERROR: cannot create the directory: res_dir" << endl;
		return 1;
	}
	string prefixName = "res_dir/res.dist.";
	
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
							double distance1 = vkssd[i]->distance(vkssd[j]);
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



