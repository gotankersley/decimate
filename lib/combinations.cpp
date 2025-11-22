#include "combinations.h"


uint64_t comb(int n, int k) { //Binomial co-efficient of n-choose-k
	//IMPORTANT NOTE: This only handles values in the range of uint64_t	
	//Larger values can use: void fmpz_bin_uiui(fmpz_t f, ulong n, ulong k)
	uint64_t res = 1;

	// Since C(n, k) = C(n, n-k), we can optimize by choosing the smaller k
	if (k > n - k) {
		k = n - k;
	}

    // Calculate value of [n * (n-1) * ... * (n-k+1)] / [k * (k-1) * ... * 1]
    for (int i = 0; i < k; ++i) {
        res *= (n - i);
        res /= (i + 1);
    }
    return res;
}

uint64_t comb_rank(std::vector<uint8_t> vals) {
	int k = vals.size();
	uint64_t rank = 0;
	
	for (int i = 0; i < k; i++) {
		rank += comb(vals[k-1-i], k-i);
	}
	return rank;
}

std::vector<uint8_t> comb_unrank(uint64_t rank, int n, int k) {
	std::vector<uint8_t> vals(k);
	
	for (int i = 0; i < k; i++) {
		while (comb(n, k-i) > rank) {
			n--;
		}
		vals[k-i-1] = n;
		rank -= comb(n, k - i);
	}
	return vals;
}