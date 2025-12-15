#include <iostream>
#include <cstdint>
#include <vector>
#include <cassert>

#include "flint/fmpz.h"
#include "flint/arith.h"
#include "../lib/io_lib.h"
#include "../lib/combinations.h"

using std::cout, std::endl;


const int N = 10;
const int K = 3;


int main() {
	int combs = comb(N, K);
	//Small
	cout << "Testing Small functions:" << endl;
	for (int r = 0; r < combs; r++) {
		std::vector<uint8_t> combVals(K);		
		comb_unrank(r, N, K, combVals);
		uint64_t rank = comb_rank(combVals);
		assert(r == rank && "rank does not match!");
		cout << "Unranking: " << r << endl;
		printVector(combVals);
	}
	
	//Big
	fmpz_t bigR;
	fmpz_t bigRank;
	
	fmpz_init(bigR);
	fmpz_init(bigRank);
	cout << "Testing Big functions (but with small values):" << endl;
	for (int r = 0; r < combs; r++) {
		fmpz_set_ui(bigR, r);
		std::vector<uint8_t> combVals(K);		
		comb_unrank(bigR, N, K, combVals);
		
		comb_rank(combVals, bigRank);
		int rank = fmpz_get_ui(bigRank);
		assert(r == rank && "rank does not match!");
		cout << "Unranking: " << r << endl;		
		printVector(combVals);
	}
	fmpz_clear(bigR);
	fmpz_clear(bigRank);
	
	return 0;
}