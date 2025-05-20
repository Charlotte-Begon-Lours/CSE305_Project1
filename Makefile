CXX = g++
CFLAGS = -pthread -std=c++14 -Wall
NVCC = /usr/local/cuda/bin/nvcc

testing: testing.cpp
	$(CXX) $(CFLAGS) -o testing testing.cpp

hello: hello.cpp
	$(CXX) $(CFLAGS) -o hello hello.cpp

clean:
	rm -f testing
	rm -f hello