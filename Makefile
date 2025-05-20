CXX = g++
CFLAGS = -pthread -std=c++14 -Wall
NVCC = /usr/local/cuda/bin/nvcc

testing: testing.cpp
	$(CXX) $(CFLAGS) -o testing testing.cpp

main: main.cpp
	$(CXX) $(CFLAGS) -o main main.cpp

clean:
	rm -f testing
	rm -f main