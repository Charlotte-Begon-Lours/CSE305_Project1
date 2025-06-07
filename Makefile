CXX = g++
CFLAGS = -pthread -std=c++14 -Wall
NVCC = /usr/local/cuda/bin/nvcc


simple_approach_test: simple_approach_test.cpp
	$(CXX) $(CFLAGS) -o simple_approach_test simple_approach_test.cpp

clean:
	rm -f simple_approach_test
