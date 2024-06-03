#TARGET=minhash
TARGET=main

CXX=g++
CXXFLAGS=-g -O2 -std=c++14 -fopenmp

INCLUDE=../build/include
LIBS=../build/lib
#OBJECTS=test.o 
OBJECTS=main.o 

$(TARGET):$(OBJECTS)
	$(CXX) $(CXXFLAGS)  $^ -o $@ -I$(INCLUDE) -lz -lrabbitsketch_static -L$(LIBS)
#	$(CXX) $(CXXFLAGS)  $^ -o $@ -lgsl -lgslcblas -lm

#%.o:%.cpp
#test.o:test.cpp
main.o:test_KSSD.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ -I$(INCLUDE) -lz



clean:
	rm *.o $(TARGET)
