#include <iostream>
#include <cstdint>
#include <vector>
#include <cassert>

#include "flint/fmpz.h"
#include "flint/arith.h"
#include "../lib/io_lib.h"
#include "../lib/set_partitions.h"

using std::cout, std::endl;


const int N = 6;
const int K = 3;
const int M = 2;




int main() {
	
	gen_coeff_table(N, K, M);
	cout << "Generating - N: " << N << ", K: " << K << ", M: " << M << endl;
	return 0;
}