#include <iostream>
#include <cstdint>
#include <vector>


#include "flint/fmpz.h"
#include "flint/arith.h"
#include "lib/near_entropic.h"

using std::cout, std::endl;



const int MAX_SYM = 16;
const int SEQ_LEN = 16;



int main() {
		
	std::vector<uint8_t> valSeq = {1,7,7,1,14,7,0,11,2,13,13,12,11,2,0,7};
	fmpz_t entRank;
	fmpz_init(entRank);
	near_entropic_rank(valSeq, MAX_SYM, entRank);			
	cout << "---" << endl;
	std::vector<uint8_t> rgfSeq(SEQ_LEN);
	near_entropic_unrank(entRank, SEQ_LEN, MAX_SYM, rgfSeq);	
	printVector(rgfSeq);
	fmpz_clear(entRank);
	return 0;
}