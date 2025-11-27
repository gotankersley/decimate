#include <iostream>
#include <iomanip>
#include <cstdint>
#include <vector>
#include <cstdlib>
#include <string>
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"
#include "flint/arith.h"
#include <fstream>
#include "lib/near_entropic.h"
using std::cout, std::endl;



const int SEQ_LEN = 1000;
const int MAX_SYM = 256;
const int K_SHAPING = 1;

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

	//printVector(valSeq);
	cout << "Entropy of Seq: " << measureEntropy(valSeq, MAX_SYM) << endl;
	
	fmpz_t normRank;
	fmpz_init(normRank);
	b2n(valSeq, MAX_SYM, normRank);
	//fmpz_print(normRank);
		
	std::vector<uint8_t> rgfSeq(SEQ_LEN+K_SHAPING);
	near_entropic_unrank(normRank, SEQ_LEN+K_SHAPING, MAX_SYM, rgfSeq);	
	fmpz_clear(normRank);
	
	//printVector(rgfSeq);
	cout << "Entropy of Transform: " << measureEntropy(rgfSeq, MAX_SYM) << endl;
	
	return 0;
}