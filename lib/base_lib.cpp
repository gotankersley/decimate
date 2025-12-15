#include "base_lib.h"


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


// Compute the information content
double measureEntropy(std::vector<int>& counts, int seqLen) {
	double n = seqLen;
    if (n == 0.0) return 0.0;
	
	double total = 0.0;
	for (int i = 0; i < (int)counts.size(); i++) {
		int count = counts[i];
		if (count == 0) continue;		
		total += (count * std::log2(count/n));
	}
	return -total;
	
}