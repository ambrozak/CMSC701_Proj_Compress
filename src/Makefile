CXX      = g++
CXXFLAGS = -std=c++17 -O3 -march=native -Wall -Wextra

.PHONY: all clean

all: build_index query_index build_plain_index query_plain_index 

build_index: build_index.cpp
	$(CXX) $(CXXFLAGS) -o build_index build_index.cpp

query_index: query_index.cpp
	$(CXX) $(CXXFLAGS) -o query_index query_index.cpp

build_plain_index: build_plain_index.cpp
	$(CXX) $(CXXFLAGS) -o build_plain_index build_plain_index.cpp

query_plain_index: query_plain_index.cpp
	$(CXX) $(CXXFLAGS) -o query_plain_index query_plain_index.cpp

clean:
	rm -f build_index query_index build_plain_index query_plain_index
