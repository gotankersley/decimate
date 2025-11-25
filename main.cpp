#include <iostream>
#include <cstdint>
#include <vector>


#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"
#include "flint/arith.h"
#include "lib/near_entropic.h"

using std::cout, std::endl;



const int MAX_SYM = 16;
const int SEQ_LEN = 16;




int main() {
		
	std::vector<uint8_t> valSeq = {1,7,7,1,14,7,0,11,2,13,13,12,11,2,0,7};
	//double ent = measureEntropy(valSeq, MAX_SYM);
	//cout << "ent: " << ent << endl;
	/*
	fmpz_mat_t table;
	gen_rgf_table(16, 8, table);
	serialize_mat("8.mat", table);
	fmpz_mat_clear(table);
	
	fmpz_mat_t table2;
	deserialize_mat("8.mat", table2);
	fmpz_mat_clear(table2);
	*/
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