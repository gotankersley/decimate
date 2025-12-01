#include <iostream>
//#include <iomanip>
#include <cstdint>
#include <vector>
#include <cstdlib>
#include <string>
#include <fstream>
#include <cassert>
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"
#include "flint/arith.h"
#include "lib/near_entropic.h"
using std::cout, std::endl;



const int SEQ_LEN = 8;//10000;
const int MAX_SYM = 8;
const int K_SHAPING = 4;

void b2n(std::vector<uint8_t>& digits, int base, fmpz_t nOut) {
	fmpz_zero(nOut);
	for (int i = 0; i < (int)digits.size(); i++) {
		uint8_t digit = digits[i];
		fmpz_t power;
		fmpz_init(power);
		fmpz_ui_pow_ui(power, base, i);
		fmpz_addmul_ui(nOut, power, digit);
		fmpz_clear(power);
	}
	
}

int main() {
	srand(42); 	
	std::vector<uint8_t> valSeq(SEQ_LEN);
	
	
	for (int i = 0; i < SEQ_LEN; i++) {
		valSeq[i] = rand()%MAX_SYM;
	}
	//cout << "Val Seq: " << endl;
	//printVector(valSeq);
	//cout << "Entropy of Seq: "<< std::fixed << std::setprecision(2) << measureEntropy(valSeq, MAX_SYM) << endl;
	
	fmpz_t normRank;
	fmpz_init(normRank);
	b2n(valSeq, MAX_SYM, normRank);

		
	fmpz_t normCopy;
	fmpz_init(normCopy);
	fmpz_set(normCopy, normRank);
	std::vector<uint8_t> rgfSeq(SEQ_LEN+K_SHAPING);
	near_entropic_unrank(normRank, SEQ_LEN+K_SHAPING, MAX_SYM, rgfSeq);		
	fmpz_clear(normRank);
	
	//cout << "RGF Seq: " << endl;
	//printVector(rgfSeq);
	//cout << "Entropy of Transform: " << std::fixed << std::setprecision(2) << measureEntropy(rgfSeq, MAX_SYM) << endl;
	
	fmpz_t rank;	
	fmpz_init(rank);
	near_entropic_rank(rgfSeq, MAX_SYM, rank);				
	//	
	assert(fmpz_equal(rank, normCopy) && "rank does not match!");				
	fmpz_clear(rank);
	fmpz_clear(normCopy);
	
	return 0;
}