CXX = g++
CFLAGS = -pthread -std=c++14 -Wall
NVCC = /usr/local/cuda/bin/nvcc

testing: testing.cpp
	$(CXX) $(CFLAGS) -o testing testing.cpp

main: main.cpp
	$(CXX) $(CFLAGS) -o main main.cpp

seq_fanny: seq_fanny.cpp
	$(CXX) $(CFLAGS) -o seq_fanny seq_fanny.cpp

clean:
	rm -f testing
	rm -f main
	rm -f seq_fanny