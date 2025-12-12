#include <iostream>
#include <cstdint>
#include <vector>
#include <cassert>

#include "flint/fmpz.h"
#include "flint/arith.h"
#include "../lib/rgf.h"

using std::cout, std::endl;


const int N = 5;
const int K =3;


void printVector(std::vector<uint8_t>& vals);// {    
//	for (int val : vals) {
//		cout << +val << ","; 
//	}
//	cout << endl;
//}


int main() {
		
	fmpz_t bigStir;
	fmpz_init(bigStir);
	arith_stirling_number_2(bigStir, N, K);
	uint64_t stir2 = fmpz_get_ui(bigStir);
	fmpz_clear(bigStir);
	
	for (uint64_t r = 0; r < stir2; r++) {	
		fmpz_t stirRank;
		fmpz_init(stirRank);
		fmpz_set_ui(stirRank, r);
		
		std::vector<uint8_t> rgfSeq(N);
		rgf_unrank(stirRank, N, K, rgfSeq);	
		fmpz_clear(stirRank);
		
		cout << "RGF: ";
		printVector(rgfSeq); 
				
		fmpz_t rank;	
		fmpz_init(rank);
		rgf_rank(rgfSeq, K, rank);				
		
		assert(fmpz_equal_ui(rank, r) && "rank does not match!");				
		fmpz_clear(rank);
	}
		
	
	return 0;
}