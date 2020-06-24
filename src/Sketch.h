#ifndef Sketch_h
#define Sketch_h

#include <map>
#include <vector>
#include <string>
//#include <string.h>
#include <stdint.h>

#include "MinHash.h"
#include "histoSketch.h"
#include "HyperLogLog.h"


/// \brief Sketch namespace
namespace Sketch{

	typedef uint64_t hash_t;

	struct Reference
	{
		// no sequence for now
		std::string name;
		std::string comment;
		uint64_t length;
		HashList hashesSorted;
		std::vector<uint32_t> counts;
	};

/// Sketching seqeunces using minhash method
	class MinHash
	{

		public:
			/// minhash init with parameters
			MinHash();

			/// minhash is updatable with multiple sequences
			void update(char * seq);

			/// merge two minhashes
			void merge(MinHash& msh);

			/// return the jaccard index
			double jaccard(MinHash * msh);			

			/// return distance defined in Mash and RabbitMash
			double distance(MinHash * msh);

			//HashList & getHashList(); //TODO: return to vector instead of HashList
			void getMinHash(); //return hash values to a vector?

			void printMinHashes();

			/// return totalSeqence length, including multiple updates
			uint64_t getTotalLength(){return totalLength;}

			// parameters

			/// set kmer size for sketching
			void setKmerSize(int kmerSizeNew) { kmerSize = kmerSizeNew; }

			/// set alphabet size: default 4 for nucleotide sequences
			void setAlphabetSize(int alphabetSizeNew) { alphabetSize = alphabetSizeNew; }

			/// set whether to preserve upper or lower case 
			void setPreserveCase(bool caseNew) { preserveCase = caseNew; }

			/// set whether to user 64bit hashes: default true
			/// use 32 bit hashes is use64 is false
			void setUse64(bool use64New) { use64 = use64New; }

			/// set hash seed
			void setSeed(uint32_t seedNew) { seed = seedNew; }

			/// set sketch size: default 1000
			void setSketchSize(uint64_t sketchSizeNew) { sketchSize = sketchSizeNew; }

			/// set whether to use reverse complement for each kmer
			void setNoncanonical(bool noncanonicalNew) { noncanonical = noncanonicalNew; }

			/// return kmerSize
			int getKmerSize() { return kmerSize; }

			/// return alphabet size
			uint32_t getAlphabetSize() { return alphabetSize; }

			/// return whether to preserve case
			bool isPreserveCase() { return preserveCase; }

			/// return whether to use 64bit hash
			bool isUse64() { return use64; }

			/// return hash seed
			uint32_t getSeed() {return seed; }

			/// return sketch size
			uint32_t getSketchSize() {return sketchSize; }

			/// return whether to use reverse complement
			bool isNoncanonical() { return noncanonical; }

		private:
			bool needToList = true;
			double kmerSpace;
			//Parameters parameters;
			MinHashHeap * minHashHeap;
			Reference reference;
			uint64_t totalLength = 0;
			/// get pValue for distance
			double pValue(uint64_t x, uint64_t lengthRef, uint64_t lengthQuery, double kmerSpace, uint64_t sketchSize);
			void heapToList();

			//parameters
			int kmerSize = 21;
			uint32_t alphabetSize = 4; //nuc sequences
			bool preserveCase = false;
			bool use64 = true;
			uint32_t seed = 42;
			uint64_t sketchSize = 1000; //minHashesPerWindow
			bool noncanonical = false;
	};

	class WMinHash{
		public:
			//WMinHash(Parameters parametersNew);
			WMinHash();
			~WMinHash();

			void update(char * seq);
			double wJaccard(WMinHash * wmh);
			double distance(WMinHash * wmh);
			void getWMinHash();

			void setKmerSize(int kmerSizeNew) { kmerSize = kmerSizeNew; }
			void setAlphabetSize(int alphabetSizeNew) { alphabetSize = alphabetSizeNew; }
			void setNumBins(int numBinsNew) { numBins = numBinsNew; }
			void setMinimizerWindowSize(int minimizerWindowSizeNew) { minimizerWindowSize = minimizerWindowSizeNew; }
			void setHistoSketchSize(int histoSketchSizeNew); //{ histoSketchSize = histoSketchSizeNew; }
			void setHistoDimension(int histoDimensionNew); //{ histoDimension = histoDimensionNew; }
			void setParaDecayWeight(double paraDecayWeightNew) { paraDecayWeight = paraDecayWeightNew; }
			void setApplyConceptDrift(bool applyConceptDriftNew) { applyConceptDrift = applyConceptDriftNew; }

			int getKmerSize() { return kmerSize; }
			int getAlphabetSize() { return alphabetSize; }
			int getNumBins() { return numBins; }
			int getMinimizerWindowSize() { return minimizerWindowSize; }
			int getHistoSketchSize() { return histoSketchSize; }
			int getHistoDimension() { return histoDimension; }
			double getParaDecayWeight() { return paraDecayWeight; }
			bool isApplyComceptDrift() { return applyConceptDrift; }


			

		private:
			//Parameters parameters;
			int kmerSize = 21;
			int alphabetSize = 4;
			int numBins = 194481;
			int minimizerWindowSize = 9;
			int histoSketchSize = 50;
			int histoDimension = 194481;
			double paraDecayWeight = 0.0;
			bool applyConceptDrift = false;
			double decayWeight = 1.0;

			bool needToCompute = true;
			double * binsArr;
			double * countMinSketch; 
			std::vector<uint64_t> sketches;
			std::vector<Bin> kmerSpectrums;
			uint32_t * histoSketches;
			double * histoWeight;

			void computeHistoSketch();

			double * r;
			double * c;
			double * b;



	};

	//OMinHash
	struct OSketch {
		//std::string       name;
		int               k, l, m;
		std::vector<char> data;
		std::vector<char> rcdata;

		bool operator==(const OSketch& rhs) const {
			return k == rhs.k && l == rhs.l && m == rhs.m && data == rhs.data && rcdata == rhs.rcdata;
		}

	};

	/// Sketching and compare sequences or strings using Order MinHash algorithm.
	class OMinHash{

		public:
			/// Order MinHash constructor
			OMinHash(){};
			/// Order MinHash constructor for sketching sequences using default parameters
			OMinHash(char * seqNew);
			~OMinHash() {if (rcseq != NULL) delete rcseq;};

			/// return sketch result in OSketch type
			OSketch getSektch(){ return sk;}

			/// build a order minhash Sketch
			/// seqNew is NULL pointer in default
			/// if seqNew is NULL pointer, buildSketch() will rebuild sketh using old data.
			/// This is useful when chaning parameters and build a new sketch
			void buildSketch(char * seqNew);

			/// return similarity (not jaccard index) between two Order MinHash sketches 
			double similarity(OMinHash & omh2);

			/// return distance between tow Order MinHash sketches
			double distance(OMinHash & omh2)
			{
				return (double)1.0 - similarity(omh2);
			}

			/// set kmer size: default 21
			void setK(int k){ m_k = k; }

			/// set "l": default 2 (normally 2 - 5)
			void setL(int l){ m_l = l; }

			/// set "m": default 500
			void setM(int m){ m_m = m; }

			/// set seed value: default 32
			void setSeed(uint64_t seedNew) { mtSeed = seedNew; }

			/// choose whether to deal with reverse complement sequences: default false 
			/// reverse complement is normally used in biological sequences such as DNA or protein
			void setReverseComplement(bool isRC){rc = isRC;}

			/// return kmer size
			int getK(){return m_k;}	

			/// return "l" value
			int getL(){return m_l;}	

			/// return "m" value
			int getM(){return m_m;}		

			/// return seed value
			uint64_t getSeed() { return mtSeed; }

			/// return whether to deal with reverse complement sequences
			bool isReverseComplement(){return rc;}

		private:

			char * seq = NULL;
			char * rcseq = NULL;
			//Parameters parameters;
			int m_k = 21, m_l = 2, m_m = 500;
			//reverse complement
			bool rc = false;
			OSketch sk;
			uint64_t mtSeed = 32; //default value

			void sketch();

			inline void compute_sketch(char * ptr, const char * seq);

			double compare_sketches(const OSketch& sk1, const OSketch& sk2, 
											  ssize_t m = -1, bool circular = false);
			double compare_sketch_pair(const char* p1, const char* p2,
									   unsigned m, unsigned k, unsigned l, bool circular);

	};
	
	class HyperLogLog{
		
		public:
			HyperLogLog(int np):core_(1uL<<np,0),np_(np),is_calculated_(0),estim_(EstimationMethod::ERTL_MLE),jestim_(JointEstimationMethod::ERTL_JOINT_MLE) {};
			~HyperLogLog(){};
			void update(const char* seq);
			HyperLogLog merge(const HyperLogLog &other) const;
			void printSketch();
			double distance(const HyperLogLog &h2) const {return jaccard_index(h2);}
			double jaccard_index(HyperLogLog &h2); 
			double jaccard_index(const HyperLogLog &h2) const; 

		protected:
			std::vector<uint8_t> core_;//sketchInfo; 
			mutable double value_; //cardinality
			uint32_t np_; // 10-20
			mutable uint8_t is_calculated_;
			EstimationMethod                        estim_;
			JointEstimationMethod                  jestim_;
			//HashStruct                                 hf_;
		
		private:
			uint32_t p() const {return np_;}//verification
			uint32_t q() const {return (sizeof(uint64_t) * CHAR_BIT) - np_;}
			uint64_t m() const {return static_cast<uint64_t>(1) << np_;}
			size_t size() const {return size_t(m());}
			bool get_is_ready() const {return is_calculated_;}
			const auto &core()    const {return core_;}
			EstimationMethod get_estim()       const {return  estim_;}
			JointEstimationMethod get_jestim() const {return jestim_;}
			void set_estim(EstimationMethod val) { estim_ = std::max(val, ERTL_MLE);}
			void set_jestim(JointEstimationMethod val) { jestim_ = val;}
			void set_jestim(uint16_t val) {set_jestim(static_cast<JointEstimationMethod>(val));}
			void set_estim(uint16_t val)  {estim_  = static_cast<EstimationMethod>(val);}

			// Returns cardinality estimate. Sums if not calculated yet.
			double creport() const {
				csum();
				return value_;
			}
			double report() noexcept {
				csum();
				return creport();
			}


			//private:
			void add(uint64_t hashval);
			void addh(const std::string &element);
			double alpha()          const {return make_alpha(m());}
			static double small_range_correction_threshold(uint64_t m) {return 2.5 * m;}
			double union_size(const HyperLogLog &other) const;
			// Call sum to recalculate if you have changed contents.
			void csum() const { if(!is_calculated_) sum(); }
			void sum() const {
				const auto counts(sum_counts(core_)); // std::array<uint32_t, 64>  // K->C
				value_ = calculate_estimate(counts, estim_, m(), np_, alpha(), 1e-2);
				is_calculated_ = 1;
			}
			std::array<uint32_t,64> sum_counts(const std::vector<uint8_t> &sketchInfo) const;
			double calculate_estimate(const std::array<uint32_t,64> &counts, EstimationMethod estim, uint64_t m, uint32_t p, double alpha, double relerr) const; 
			template<typename T>
				void compTwoSketch(const std::vector<uint8_t> &sketch1, const std::vector<uint8_t> &sketch2, T &c1, T &c2, T &cu, T &cg1, T &cg2, T &ceq) const;
			template<typename T>
				double ertl_ml_estimate(const T& c, unsigned p, unsigned q, double relerr=1e-2) const; 
			template<typename HllType>
				std::array<double, 3> ertl_joint(const HllType &h1, const HllType &h2) const; 



	};


}//namespace sketch

#endif //Sketch_h
