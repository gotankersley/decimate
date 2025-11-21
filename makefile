CXX = g++
CXXFLAGS = -Wall -std=c++17 
LDFLAGS = -lflint # Library linking

build: 
	$(CXX) $(CXXFLAGS) -o decimate main.cpp $(LDFLAGS)


clean:
	rm -f *.o decimate