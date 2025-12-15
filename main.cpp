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
#include "lib/base_lib.h"
#include "lib/near_entropic.h"
using std::cout, std::endl;

const int SEQ_LEN = 100;
const int MAX_SYM = 16;
const int K_SHAPING = 1;
const int K_SEARCH = 10;
const bool VERIFY = false;



int main() {
	srand(44); 	
	std::vector<uint8_t> valSeq(SEQ_LEN);
	
	//Generate a random sequence
	for (int i = 0; i < SEQ_LEN; i++) {
		valSeq[i] = rand()%MAX_SYM;
	}
	
	
	fmpz_t normRank;
	fmpz_init(normRank);
	std::vector<int> valCounts(MAX_SYM);
	b2n(valSeq, MAX_SYM, valCounts, normRank);	
	double valEntropy = measureEntropy(valCounts, valSeq.size());
	cout << "Entropy of Val Seq: "<< std::fixed << std::setprecision(2) << valEntropy << endl;
	
	fmpz_t normCopy;
	fmpz_init(normCopy);
	int bestK = -1;
	double bestDelta = -1;
	for (int k = 0; k < K_SEARCH; k++) {
		fmpz_set(normCopy, normRank);
		std::vector<uint8_t> rgfSeq(SEQ_LEN+k);
		std::vector<int> entCounts(MAX_SYM);
		near_entropic_unrank(normCopy, SEQ_LEN+k, MAX_SYM, entCounts, rgfSeq);		
		double entEntropy = measureEntropy(entCounts, rgfSeq.size());
		double delta = valEntropy - entEntropy;
		if (delta > bestDelta) {
			bestDelta = delta;
			bestK = k;
		}
		cout << " - Entropy of Ent Seq " << k << ": " << std::fixed << std::setprecision(2) << entEntropy << " | Delta: " << delta << endl;	
	}
	cout << "Best K: " << bestK << ", Delta: " << bestDelta << endl;
	fmpz_clear(normCopy);
			
	//if (VERIFY) {
	//	fmpz_t rank;	
	//	fmpz_init(rank);
	//	near_entropic_rank(rgfSeq, MAX_SYM, rank);
	//
	//	assert(fmpz_equal(rank, normRank) && "rank does not match!");				
	//	fmpz_clear(rank);
	//}
	fmpz_clear(normRank);
	
	return 0;
}