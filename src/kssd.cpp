#include "Sketch.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <math.h>
#include <string.h>


namespace Sketch
{
	const unsigned int primer[25] = 
	{
		251, 509, 1021, 2039, 4093, 8191, 16381,
		32749, 65521, 131071, 262139, 524287,
		1048573, 2097143, 4194301, 8388593, 16777213,
		33554393, 67108859, 134217689, 268435399,
		536870909, 1073741789, 2147483647, 4294967291
	};
	#define H1(K, HASH_SZ) ((K) % (HASH_SZ))
	#define H2(K, HASH_SZ) (1 + (K) % ((HASH_SZ)-1))
	#define HASH(K, I, HASH_SZ) ((H1(K, HASH_SZ) + I * H2(K, HASH_SZ)) % HASH_SZ)

	int * KSSDParameters::shuffleN(int n, int base)
	{
		int * arr;
		arr = (int* ) malloc(n * sizeof(int));
		for(int i = 0; i < n; i++){
			arr[i] = i + base;
		}
		
		return shuffle(arr, n);
	}
	
	int * KSSDParameters::shuffle(int arr[], int length)
	{
		if(length > RAND_MAX){
			fprintf(stderr, "shuffling array length %d must be less than RAND_MAX: %d", length, RAND_MAX);
			exit(1);
		}
		//srand(time(NULL));
		srand(23);
		int j, tmp;
		for(int i = length-1; i > 0; i--)
		{
			j = rand() % (i + 1);
			tmp = arr[i];
			arr[i] = arr[j];
			arr[j] = tmp;
		}
		
		return arr;
	}
	
	int KSSDParameters::get_hashSize(int half_k, int drlevel)
	{
		int dim_reduce_rate = 1 << 4 * drlevel;
		uint64_t ctx_space_sz = 1LLU << 4 * (half_k - drlevel);
		int primer_id = 4 * (half_k - drlevel) - CTX_SPC_USE_L - 7;
		if(primer_id < 0 || primer_id > 24)
		{
			int k_add = primer_id < 0 ? (1 + (0 - primer_id) / 4) : -1 * (1 + (primer_id - 24) / 4);
			fprintf(stderr, "get_hashSize(): primer_id: %d out of range(0 ~ 24), by formula:\n"
									"int primer_id = 4 * (half_k - drlevel) - CTX_SPC_USE_L - 7;\n"
									"this might caused by too small or too large k\n"
									"half kmer length = %d\n"
									"dim reduction level %d\n"
									"ctx_space size = %llu\n"
									"try rerun the program with option -k = %d\n",
									primer_id, half_k, drlevel, ctx_space_sz, half_k + k_add);
	
		}
		int hashSize = primer[primer_id];
		fprintf(stderr, "dimension reduced %d\n"
										"ctx_space size  %llu\n"
										"half_k is: %d\n"
										"drlevel is: %d\n"
										"primer_id is: %d\n"
										"hashSize is: %u \n",
										dim_reduce_rate, ctx_space_sz, half_k, drlevel, primer_id, hashSize);
		
		return hashSize;
	}

	int KSSD::get_half_k(){
		return half_k;
	}
	
	int KSSD::get_half_subk(){
		return half_subk;
	}
	int KSSD::get_drlevel(){
		return drlevel;
	}

	vector<uint64_t> KSSD::storeHashes(){
		return hashList;
	}

	void KSSD::loadHashes(vector<uint64_t> hashArr){
		hashList = hashArr;
	}

	void KSSD::SetToList()
	{
		if(hashList.empty())//first update the hashList
		{
			for(auto x : hashSet)
				hashList.push_back(x);
		}
		else//not the first update the hashList, need merge.
		{
			for(int i = 0; i < hashList.size(); i++)
				hashSet.insert(hashList[i]);
			hashList.clear();
			for(auto x : hashSet)
			{
				hashList.push_back(x);
				//cerr << x << endl;
			}
			//cerr << endl;
		}
	
		//cerr << "hashList.size() is: " << hashList.size() << endl;
		unordered_set<uint64_t>().swap(hashSet);//release the memory of hashSet.
		std::sort(hashList.begin(), hashList.end());
	}
	
	
	void KSSD::update(char * seq)
	{
		uint64_t tuple = 0LLU, rvs_tuple = 0LLU, uni_tuple, dr_tuple, pfilter;
		int keyCount = 0;
		uint64_t base = 1;
		int length = strlen(seq);
		//cout << "---------------------length is: " << length << endl;
		for(int i = 0; i < length; i++)
		{
			char ch = seq[i];
			int basenum = BaseMap[(int)ch];
			if(basenum != DEFAULT_CHAR_KSSD)
			{
				tuple = ((tuple << 2) | basenum) & tupmask;
				rvs_tuple = (rvs_tuple >> 2) + (((uint64_t)basenum ^3LLU) << rev_addmove); 
				base++;
			}
			if(i >= kmer_size-1)
			{
				uni_tuple = tuple < rvs_tuple ? tuple : rvs_tuple;
				int dim_id = (uni_tuple & domask) >> (half_outctx_len * 2);
				pfilter = shuffled_dim[dim_id];
				
				if(pfilter >= dim_end || pfilter < dim_start) continue;
				pfilter -= dim_start;
				dr_tuple = (((uni_tuple & undomask0) | ((uni_tuple & undomask1) << (kmer_size * 2 - half_outctx_len * 4))) >> (drlevel * 4)) | pfilter; 
				hashSet.insert(dr_tuple);
			}//end if i > kmer_size
		}//end for, of a sequence
	
		SetToList();
		//cerr << "the sketch size of kssd is: " << hashList.size() << endl;
		if(hashList.size() > hashLimit)
			fprintf(stderr, "the context space is too crowd, try rerun the program using -k %d\n", half_k + 1);
	}//end KSSD::update
	
	double KSSD::jaccard(KSSD* kssd)
	{
		double jaccard = 1.0;
		int size1 = hashList.size();
		int size2 = kssd->hashList.size();
		int i = 0, j = 0;
		int common = 0;
		while(i < size1 && j < size2)
		{
			if(hashList[i] < kssd->hashList[j])
				i++;
			else if(hashList[i] > kssd->hashList[j])
				j++;
			else{
				i++;
				j++;
				common++;
			}
		}
		int denom = size1 + size2 -common;
		jaccard = (double)common / (double)denom;
		return jaccard;
	}
	
	double KSSD::distance(KSSD* kssd)
	{
		double distance = 1.0;
		double jaccard_ = jaccard(kssd);
		if(jaccard_ == 0.0)
			distance = 1.0;
		else if (jaccard_ == 1.0)
			distance = 0.0;
		else
			distance = -1.0 / (double)kmer_size * log((2 * jaccard_) / (1.0 + jaccard_));
		return distance;
	}
	
	void KSSD::printHashes()
	{
		for(int i = 0; i < hashList.size(); i++)
			cout << hashList[i] << endl;
	}



}
