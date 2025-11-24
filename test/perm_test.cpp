#include <iostream>
#include <cstdint>
#include <vector>
#include <cassert>

#include "flint/fmpz.h"
#include "flint/arith.h"
#include "../lib/permutations.h"

using std::cout, std::endl;


const int P = 3;


void printVector(std::vector<uint8_t>& vals) {    
	for (int val : vals) {
		cout << val << ","; 
	}
	cout << endl;
}

int main() {
	
	fmpz_t bigFact;
	fmpz_fac_ui(bigFact, P);
	uint64_t fact = fmpz_get_ui(bigFact);
	
	for (int r = 0; r < fact; r++) {
		fmpz_t rank;	

		std::vector<uint8_t> perm(P);
		myrvold_unrank(perm, rank);
		cout << "Unrank: ";
		printVector(uperm); 
		
		std::vector<uint8_t> perm = {2, 0, 1};
		myrvold_rank(mrank, perm);
		cout << "Perm: ";
		printVector(perm);
		fmpz_print(mrank);
		cout << endl;
		
		assert(r == rank && "rank does not match!");
	}
	int combs = comb(N, K);
	
	
	return 0;
}