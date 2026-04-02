CXX = g++
CXXFLAGS = -O3 -std=c++17 -march=native -ffast-math
LDFLAGS = -lm

TARGETS = compress query

all: $(TARGETS)

compress: compress.cpp
	$(CXX) $(CXXFLAGS) -o compress compress.cpp $(LDFLAGS)

query: query.cpp
	$(CXX) $(CXXFLAGS) -o query query.cpp $(LDFLAGS)

clean:
	rm -f $(TARGETS)

.PHONY: all clean