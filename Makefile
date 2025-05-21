CXX = g++
CFLAGS = -pthread -std=c++14 -Wall
NVCC = /usr/local/cuda/bin/nvcc

testing: testing.cpp
	$(CXX) $(CFLAGS) -o testing testing.cpp

main: main.cpp
	$(CXX) $(CFLAGS) -o main main.cpp

seq_fanny: seq_fanny.cpp
	$(CXX) $(CFLAGS) -o seq_fanny seq_fanny.cpp

simple_approach_test: simple_approach_test.cpp
	$(CXX) $(CFLAGS) -o simple_approach_test simple_approach_test.cpp

clean:
	rm -f testing
	rm -f main
	rm -f seq_fanny
	rm -f simple_approach_test