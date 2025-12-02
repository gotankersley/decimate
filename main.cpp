#include <iostream>
#include <iomanip>
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

const int SEQ_LEN = 10000;
const int MAX_SYM = 128;
const int K_SHAPING = 0;
const bool VERIFY = false;

void b2n(std::vector<uint8_t>& digits, int base, std::vector<int>& countsOut, fmpz_t nOut) {
	fmpz_zero(nOut);
	for (int i = 0; i < (int)digits.size(); i++) {
		uint8_t digit = digits[i];
		fmpz_t power;
		fmpz_init(power);
		fmpz_ui_pow_ui(power, base, i);
		fmpz_addmul_ui(nOut, power, digit);
		fmpz_clear(power);
		countsOut[digit]++;
	}
}

int main() {
	srand(42); 	
	std::vector<uint8_t> valSeq(SEQ_LEN);
	
	//Generate a random sequence
	for (int i = 0; i < SEQ_LEN; i++) {
		valSeq[i] = rand()%MAX_SYM;
	}
	
	
	fmpz_t normRank;
	fmpz_init(normRank);
	std::vector<int> valCounts(MAX_SYM);
	b2n(valSeq, MAX_SYM, valCounts, normRank);	
		
	fmpz_t normCopy;
	fmpz_init(normCopy);
	fmpz_set(normCopy, normRank);
	std::vector<uint8_t> rgfSeq(SEQ_LEN+K_SHAPING);
	std::vector<int> entCounts(MAX_SYM);
	near_entropic_unrank(normRank, SEQ_LEN+K_SHAPING, MAX_SYM, entCounts, rgfSeq);		
	fmpz_clear(normRank);
	
	double valEntropy = measureEntropy(valCounts, SEQ_LEN);
	double entEntropy = measureEntropy(entCounts, SEQ_LEN);
	cout << "Entropy of Val Seq: "<< std::fixed << std::setprecision(2) << valEntropy << endl;
	cout << "Entropy of Ent Seq: "<< std::fixed << std::setprecision(2) << entEntropy << endl;
	cout << "Delta: "<< std::fixed << std::setprecision(2) << valEntropy - entEntropy << endl;

	if (VERIFY) {
		fmpz_t rank;	
		fmpz_init(rank);
		near_entropic_rank(rgfSeq, MAX_SYM, rank);
	
		assert(fmpz_equal(rank, normCopy) && "rank does not match!");				
		fmpz_clear(rank);
	}
	fmpz_clear(normCopy);
	
	return 0;
}