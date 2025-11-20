CXX = g++
CXXFLAGS = -Wall -std=c++17 # Flags
#LDFLAGS = -lmylib # Library linking

build: 
	$(CXX) $(CXXFLAGS) -o decimate main.cpp


clean:
	rm -f *.o decimate