#include <iostream>
#include <cstdint>
#include <vector>
#include <cassert>
#include <cmath>
#include "flint/fmpz.h"
#include "flint/arith.h"
#include "../lib/nearer_entropic.h"

using std::cout, std::endl;


const int MAX_SYM = 4;//5;
const int SEQ_LEN = 3;

int main() {
		       
	//std::vector<uint8_t> seq = {0, 1, 2, 1, 2};				
	std::vector<uint8_t> seq = {1,0,3,3,0,1};
	fmpz_t rank;	
	fmpz_init(rank);
	nearer_entropic_rank(seq, MAX_SYM, rank);
	fmpz_clear(rank);
	/*
	uint64_t total = (uint64_t)pow(MAX_SYM, SEQ_LEN);
	
	
	for (uint64_t r = 0; r < total; r++) {		
		fmpz_t bigR;
		fmpz_init(bigR);
		fmpz_set_ui(bigR, r);	
				
		std::vector<uint8_t> seq(SEQ_LEN);
		std::vector<int> counts(MAX_SYM);
		near_entropic_unrank(bigR, SEQ_LEN, MAX_SYM, counts, seq);	
		fmpz_clear(bigR);
		
		cout << "Seq: ";
		printVector(seq); 
				
		fmpz_t entRank;	
		fmpz_init(entRank);
		near_entropic_rank(seq, MAX_SYM, entRank);
		
		assert(fmpz_equal_ui(entRank, r) && "rank does not match!");
				
		fmpz_clear(entRank);
		
	}
	*/
		
	
	return 0;
}