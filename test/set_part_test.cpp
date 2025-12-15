#include <iostream>
#include <cstdint>
#include <vector>
#include <cassert>

#include "flint/fmpz.h"
#include "flint/arith.h"
#include "../lib/set_partitions.h"

using std::cout, std::endl;


const int N = 6;
const int K = 4;
const int M = 3;


void printVector(std::vector<uint8_t>& vals);// {    
//	for (int val : vals) {
//		cout << +val << ","; 
//	}
//	cout << endl;
//}


int main() {
	
	gen_coeff_table(N, K, M);
	return 0;
}