#include <iostream>
#include <cstdint>
#include <vector>
#include <cassert>

#include "flint/fmpz.h"
#include "flint/arith.h"
#include "../lib/combinations.h"

using std::cout, std::endl;


const int N = 10;
const int K = 3;

void printVector(std::vector<uint8_t>& vals) {    
	for (int val : vals) {
		cout << val << ","; 
	}
	cout << endl;
}

int main() {
	int combs = comb(N, K);
	for (int r = 0; r < combs; r++) {
		std::vector<uint8_t> combVals(K);		
		comb_unrank(combVals, r, N, K);
		uint64_t rank = comb_rank(combVals);
		assert(r == rank && "rank does not match!");
		cout << "Unranking: " << r << endl;
		printVector(combVals);
	}
	
	
	return 0;
}