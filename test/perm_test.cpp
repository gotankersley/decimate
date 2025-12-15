#include <iostream>
#include <cstdint>
#include <vector>
#include <cassert>

#include "flint/fmpz.h"
#include "flint/arith.h"
#include "../lib/io_lib.h"
#include "../lib/permutations.h"

using std::cout, std::endl;


const int P = 5;


uint64_t factorial(int n) {
    if (n < 0) {
        // Factorial is not defined for negative numbers
        return 0; // Or throw an exception
    }
    uint64_t result = 1;
    for (int i = 1; i <= n; ++i) {
        result *= i;
    }
    return result;
}



int main() {
		
	uint64_t fact = factorial(P);
	
	for (uint64_t r = 0; r < fact; r++) {	
		fmpz_t bigR;
		fmpz_init(bigR);
		fmpz_set_ui(bigR, r);						
		std::vector<uint8_t> perm(P);
		myrvold_unrank(bigR, perm);
		fmpz_clear(bigR);
		
		cout << "Perm: ";
		printVector(perm); 
				
		fmpz_t rank;	
		fmpz_init(rank);
		myrvold_rank(perm, rank);				
		
		assert(fmpz_equal_ui(rank, r) && "rank does not match!");
		
		
		fmpz_clear(rank);
	}
		
	
	return 0;
}