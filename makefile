CXX = g++
CXXFLAGS = -Wall -std=c++17 
LDFLAGS = -lflint # Library linking

SRCS = main.cpp combinations.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = decimate

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LDFLAGS)

clean:
	rm -f *.o decimate