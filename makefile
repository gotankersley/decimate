CXX = g++ #-O3
CXXFLAGS = -Wall -std=c++17 -I./lib
LDFLAGS = -lflint #-lgmp -lgmpxx # Library linking

# --- Main Application Files ---
SRCS = main.cpp lib/near_entropic.cpp lib/nearer_entropic.cpp lib/combinations.cpp lib/permutations.cpp lib/rgf.cpp lib/set_partitions.cpp lib/io_lib.cpp lib/base_lib.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = decimate

# --- Test Files ---
# Finds all .cpp files in the ./test/ directory
TEST_SRCS = test/nearer_ent_test.cpp
#TEST_SRCS = test/set_part_test.cpp
TEST_OBJS = $(TEST_SRCS:.cpp=.o)
TEST_TARGET = run_tests


# --- Rules ---
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LDFLAGS)


# Test Link
# We use $(filter-out main.o, $(OBJS)) to link 'combinations.o' 
# but exclude 'main.o' to avoid conflict with the test runner's main().
$(TEST_TARGET): $(TEST_OBJS) $(OBJS)
	$(CXX) -o $@ $(TEST_OBJS) $(filter-out main.o, $(OBJS)) $(LDFLAGS)

# Helper target to build AND run tests
test: $(TEST_TARGET)
	./$(TEST_TARGET)


%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LDFLAGS)

clean:
	rm -f *.o decimate lib/*.o test/*.o $(TEST_TARGET)
