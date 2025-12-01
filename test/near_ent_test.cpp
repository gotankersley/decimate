#include <iostream>
#include <cstdint>
#include <vector>
#include <cassert>
#include <cmath>
#include "flint/fmpz.h"
#include "flint/arith.h"
#include "../lib/near_entropic.h"

using std::cout, std::endl;


const int MAX_SYM = 5;
const int SEQ_LEN = 3;

int main() {
		
	uint64_t total = (uint64_t)pow(MAX_SYM, SEQ_LEN);
	
	
	for (uint64_t r = 0; r < total; r++) {		
		fmpz_t bigR;
		fmpz_init(bigR);
		fmpz_set_si(bigR, r);	
				
		std::vector<uint8_t> seq(SEQ_LEN);
		near_entropic_unrank(bigR, SEQ_LEN, MAX_SYM, seq);	
		fmpz_clear(bigR);
		
		cout << "Seq: ";
		printVector(seq); 
				
		fmpz_t entRank;	
		fmpz_init(entRank);
		near_entropic_rank(seq, MAX_SYM, entRank);
		
		assert(fmpz_equal_si(entRank, r) && "rank does not match!");
				
		fmpz_clear(entRank);
		
	}
	
		
	
	return 0;
}