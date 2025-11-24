#include "combinations.h"
#include <iostream>

uint64_t comb(int n, int k) { //Binomial co-efficient of n-choose-k
	//IMPORTANT NOTE: This only handles values in the range of uint64_t	
	//Larger values can use: void fmpz_bin_uiui(fmpz_t f, ulong n, ulong k)
	uint64_t res = 1;

	if ((k < 0) || (n < k)) return 0;
	if ((2*k) > n) k = n-k;
	
	if (k > 0) {
		for(int i = 0; i <= (k-1); i++) {
			res = (res * (n-i))/(i+1);
		}
	}
	return res;
}

uint64_t comb_rank(std::vector<uint8_t>& vals) {
	int k = vals.size();
	uint64_t rankOut = 0;
	
	for (int i = 0; i < k; i++) {
		rankOut += comb(vals[k-1-i], k-i);
	}
	return rankOut;
}

void comb_unrank(uint64_t rank, int n, int k, std::vector<uint8_t>& valsOut) {	
	
	for (int i = 0; i < k; i++) {
		while (comb(n, k-i) > rank) {			
			n--;
		}
		valsOut[k-i-1] = n;
		rank -= comb(n, k - i);
	}	
}