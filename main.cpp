#include <iostream>
#include <cstdint>
#include <vector>
//#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/arith.h"
using std::cout, std::endl;

#define big fmpz_t
#define sym_t uint8_t


const bool DEBUG = true;
const int MAX_SYM = 16;
const int SEQ_LEN = 16;
const int INVALID = -1;

void near_entropic_rank( big rankOut) {
	sym_t valSeq[] = {1,7,7,1,14,7,0,11,2,13,13,12,11,2,0,7};
	if (DEBUG) cout << "Ranking: " <<  endl; 
	
	//Pre-process loop file bits -> seq	
	int valToSym[MAX_SYM]; //Put val in get sym out
	for (int i = 0; i < MAX_SYM; i++) {
		valToSym[i] = INVALID;
	}
	std::vector<sym_t> vals;
	sym_t symSeq[SEQ_LEN];
		
	int symCount = INVALID;
	//Loop read in files as i16
	for (int i = 0; i < SEQ_LEN; i++){
		sym_t val = valSeq[i];
		sym_t sym;
		if (valToSym[val] == INVALID) { //First time this symbol has been seen		
			symCount++;
			valToSym[val] = symCount;
			sym = symCount;
			vals.push_back(val);
		}
		else sym = valToSym[val];
		
		symSeq[i] = sym;	
		cout << +symSeq[i] << endl;
	}
	
	if (DEBUG) cout << "Sym Seq: " << symSeq << endl;
	//Use std set to get distinct symCount
	//Create sym seqence 
	//Create map of seq -> vals
	
	
	
	// 1. Add symbol sections
	//for (int i = 0; i < symCount; i++) {		
	//	cout << "here:" << endl;
	//	//rank += countSymbolSection(seqLen, totalSymbolsAvail, i)
	//}
		
	//if (DEBUG) count << "rank after symbol section: <<  rank << endl;
	
	//arith_stirling_number_2(rank, n, k);	
}

int main() {
		
	big rank;
	near_entropic_rank( rank);	
	//fmpz_print(rank);

	return 0;
}