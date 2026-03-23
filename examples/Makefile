INCLUDE=../build/include
LIBS=../build/lib
CXX=g++
CXXFLAGS=-g -O2 -std=c++14 -fopenmp -fPIC -march=native -I $(INCLUDE)

all: exe_test_HLL exe_test_SetSketch exe_test_Kssd exe_test_MinHash exe_test_OMH exe_generate_shuf_file exe_main exe_test_ProbMinHash exe_test_OneProbMinHash exe_eval_sketch_acc

exe_test_HLL: test_HLL.cpp
	$(CXX) $(CXXFLAGS) test_HLL.cpp -o exe_test_HLL -lz -lrabbitsketch_static -L $(LIBS)

exe_test_SetSketch: test_SetSketch.cpp
	$(CXX) $(CXXFLAGS) test_SetSketch.cpp -o exe_test_SetSketch -lz -lrabbitsketch_static -L $(LIBS)

exe_test_Kssd: test_Kssd.cpp
	$(CXX) $(CXXFLAGS) test_Kssd.cpp -o exe_test_Kssd -lz -lrabbitsketch_static -L $(LIBS) 

exe_test_MinHash: test_MinHash.cpp
	$(CXX) $(CXXFLAGS) test_MinHash.cpp  -o exe_test_MinHash -lz -lrabbitsketch_static -L $(LIBS)

exe_test_OMH: test_OMH.cpp
	$(CXX) $(CXXFLAGS) test_OMH.cpp -o exe_test_OMH -lz -lrabbitsketch_static -L $(LIBS)

exe_generate_shuf_file: generate_shuf_file.cpp
	$(CXX) $(CXXFLAGS) generate_shuf_file.cpp -o exe_generate_shuf_file -lz -lrabbitsketch_static -L $(LIBS)

exe_main: main.cpp
	$(CXX) $(CXXFLAGS) main.cpp -o exe_main -lz -lrabbitsketch_static -L $(LIBS)

exe_test_ProbMinHash: test_ProbMinHash.cpp ../src/probmh.cpp
	$(CXX) $(CXXFLAGS) -I ../src test_ProbMinHash.cpp ../src/probmh.cpp -o exe_test_ProbMinHash -lz -lrabbitsketch_static -L $(LIBS)

exe_test_OneProbMinHash: test_OneProbMinHash.cpp ../src/probmh_op.cpp
	$(CXX) $(CXXFLAGS) -I ../src test_OneProbMinHash.cpp ../src/probmh_op.cpp -o exe_test_OneProbMinHash -lz -lrabbitsketch_static -L $(LIBS)

exe_eval_sketch_acc: eval_sketch_accuracy.cpp ../src/probmh.cpp ../src/probmh_op.cpp
	$(CXX) $(CXXFLAGS) -I ../src eval_sketch_accuracy.cpp ../src/probmh.cpp ../src/probmh_op.cpp -o exe_eval_sketch_acc -lz -lrabbitsketch_static -L $(LIBS)

clean:
	rm -f *.o exe_*
