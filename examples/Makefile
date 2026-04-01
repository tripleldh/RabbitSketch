INCLUDE=../build/include
# CMake places librabbitsketch_static.a in ../build; shared lib is under ../build/lib
RABBIT_LIBDIRS=-L../build/lib -L../build
CXX=g++
#
# The static library must be rebuilt without LTO after CMake change (see ../CMakeLists.txt).
# If ../build/lib/librabbitsketch_static.a was built with -flto, conda/other g++ will fail
# at link with "bytecode stream ... lto-wrapper". Fix: re-run cmake --build in ../build.
#
# Optional: make USE_LTO=1 only if you rebuilt librabbitsketch_static.a with the same $(CXX)+-flto.
USE_LTO ?= 0
CXXFLAGS=-g1 -O3 -std=c++20 -fopenmp -march=native -DNDEBUG -funroll-loops -fno-lto -I $(INCLUDE)
ifeq ($(USE_LTO),1)
CXXFLAGS := $(filter-out -fno-lto,$(CXXFLAGS))
CXXFLAGS += -flto
endif

all: exe_test_HLL exe_test_SetSketch exe_test_Kssd exe_test_MinHash exe_test_OMH exe_generate_shuf_file exe_main exe_test_ProbMinHash exe_test_OneProbMinHash exe_test_ProbKMV exe_test_TopLProbMinHash  exe_eval_sketch_acc

exe_test_HLL: test_HLL.cpp
	$(CXX) $(CXXFLAGS) test_HLL.cpp -o exe_test_HLL $(RABBIT_LIBDIRS) -lz -lrabbitsketch_static

exe_test_SetSketch: test_SetSketch.cpp
	$(CXX) $(CXXFLAGS) test_SetSketch.cpp -o exe_test_SetSketch $(RABBIT_LIBDIRS) -lz -lrabbitsketch_static

exe_test_Kssd: test_Kssd.cpp
	$(CXX) $(CXXFLAGS) test_Kssd.cpp -o exe_test_Kssd $(RABBIT_LIBDIRS) -lz -lrabbitsketch_static 

exe_test_MinHash: test_MinHash.cpp
	$(CXX) $(CXXFLAGS) test_MinHash.cpp  -o exe_test_MinHash $(RABBIT_LIBDIRS) -lz -lrabbitsketch_static

exe_test_OMH: test_OMH.cpp
	$(CXX) $(CXXFLAGS) test_OMH.cpp -o exe_test_OMH $(RABBIT_LIBDIRS) -lz -lrabbitsketch_static

exe_generate_shuf_file: generate_shuf_file.cpp
	$(CXX) $(CXXFLAGS) generate_shuf_file.cpp -o exe_generate_shuf_file $(RABBIT_LIBDIRS) -lz -lrabbitsketch_static

exe_main: main.cpp
	$(CXX) $(CXXFLAGS) main.cpp -o exe_main $(RABBIT_LIBDIRS) -lz -lrabbitsketch_static

exe_test_ProbMinHash: test_ProbMinHash.cpp ../src/probmh.cpp
	$(CXX) $(CXXFLAGS) -I ../src test_ProbMinHash.cpp ../src/probmh.cpp -o exe_test_ProbMinHash $(RABBIT_LIBDIRS) -lz -lrabbitsketch_static

exe_test_OneProbMinHash: test_OneProbMinHash.cpp ../src/probmh_op.cpp
	$(CXX) $(CXXFLAGS) -I ../src test_OneProbMinHash.cpp ../src/probmh_op.cpp -o exe_test_OneProbMinHash $(RABBIT_LIBDIRS) -lz -lrabbitsketch_static

exe_test_ProbKMV: test_ProbKMV.cpp ../src/probkmv.cpp
	$(CXX) $(CXXFLAGS) -I ../src test_ProbKMV.cpp ../src/probkmv.cpp -o exe_test_ProbKMV $(RABBIT_LIBDIRS) -lz -lrabbitsketch_static

exe_test_TopLProbMinHash: test_TopLProbMinHash.cpp ../src/probmh.cpp
	$(CXX) $(CXXFLAGS) -I ../src test_TopLProbMinHash.cpp ../src/probmh.cpp -o exe_test_TopLProbMinHash $(RABBIT_LIBDIRS) -lz -lrabbitsketch_static

exe_eval_sketch_acc: eval_sketch_accuracy.cpp ../src/probmh.cpp ../src/probmh_op.cpp ../src/probkmv.cpp
	$(CXX) $(CXXFLAGS) -I ../src eval_sketch_accuracy.cpp ../src/probmh.cpp ../src/probmh_op.cpp ../src/probkmv.cpp -o exe_eval_sketch_acc $(RABBIT_LIBDIRS) -lz -lrabbitsketch_static

clean:
	rm -f *.o exe_*
