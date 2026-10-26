CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall

# OpenMP: Apple clang needs libomp (brew install libomp); elsewhere plain -fopenmp
ifeq ($(shell uname), Darwin)
LIBOMP = $(shell brew --prefix libomp 2>/dev/null)
OMPFLAGS = -Xclang -fopenmp -I$(LIBOMP)/include -L$(LIBOMP)/lib -lomp
CXXFLAGS += -mcpu=native
else
OMPFLAGS = -fopenmp
CXXFLAGS += -march=native
endif

closeness: src/closeness.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(OMPFLAGS)

clean:
	rm -f closeness
