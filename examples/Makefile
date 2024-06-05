INCLUDE=../build/include
LIBS=../build/lib
CXX=g++
CXXFLAGS=-g -O2 -std=c++14 -fopenmp -I $(INCLUDE)

all: exe_test_HLL exe_test_KSSD exe_test_MinHash exe_test_OMH

exe_test_HLL: test_HLL.cpp
	$(CXX) $(CXXFLAGS) test_HLL.cpp -o exe_test_HLL -lz -lrabbitsketch_static -L $(LIBS)

exe_test_KSSD: test_KSSD.cpp
	$(CXX) $(CXXFLAGS) test_KSSD.cpp -o exe_test_KSSD -lz -lrabbitsketch_static -L $(LIBS) 

exe_test_MinHash: test_MinHash.cpp
	$(CXX) $(CXXFLAGS) test_MinHash.cpp  -o exe_test_MinHash -lz -lrabbitsketch_static -L $(LIBS)

exe_test_OMH: test_OMH.cpp
	$(CXX) $(CXXFLAGS) test_OMH.cpp -o exe_test_OMH -lz -lrabbitsketch_static -L $(LIBS)


clean:
	rm -f *.o exe_*
