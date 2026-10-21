CXX      ?= clang++
CXXFLAGS ?= -std=c++17 -O3 -Wall -Wextra
SRCS      = main.cpp graph.cpp closeness.cpp null_models.cpp
OBJS      = $(SRCS:.cpp=.o)

lab3: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

%.o: %.cpp *.hpp
	$(CXX) $(CXXFLAGS) -c $<

# Build with debug info and sanitizers (slow, but catches out-of-bounds bugs)
debug: CXXFLAGS = -std=c++17 -O1 -g -Wall -Wextra -fsanitize=address,undefined
debug: clean lab3

clean:
	rm -f lab3 $(OBJS)

.PHONY: clean debug
