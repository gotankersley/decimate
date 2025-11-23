CXX = g++
CXXFLAGS = -Wall -std=c++17 -I./lib
LDFLAGS = -lflint # Library linking

SRCS = main.cpp lib/combinations.cpp lib/permutations.cpp lib/set_partitions.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = decimate

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LDFLAGS)

clean:
	rm -f *.o decimate lib/*.o