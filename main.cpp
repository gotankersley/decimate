#include <iostream>
#include <cstdint>
#include <vector>

//#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/arith.h"
#include "lib/combinations.h"
#include "lib/permutations.h"
//#include "lib/set_partitions.h"
#include <array>
using std::cout, std::endl;

#define big fmpz_t


const bool DEBUG = true;
const int MAX_SYM = 16;
const int SEQ_LEN = 16;
const int INVALID = -1;


void printSeq(const uint8_t arr[], int size) {    
	cout << "[";
    for (int i = 0; i < size; ++i) {
        cout << +arr[i] << ",";
    }
    cout << "]" << endl;
}

void printVector(std::vector<uint8_t>& vals) {    
	for (int val : vals) {
		cout << val << ","; 
	}
	cout << endl;
}



//Combinatorial function to count the ways to rank:
// - Filling up a sequence of N length (e.g. N = 3)
// - Using an alphabet of A symbols available (e.g. A = {A,B,C} = 3 )
// - Using subset K of symbols (e.g. K = 2)
void addSymbolSection(big rank, int k) {
	if (k < 0 || k > SEQ_LEN || k > MAX_SYM) return;
	
	uint64_t combinations = comb(MAX_SYM, k);
	
	big factorialK;
	fmpz_fac_ui(factorialK, k);
	
	big stir2;
	arith_stirling_number_2(stir2, SEQ_LEN, k);

	big product;
	fmpz_mul(product, factorialK, stir2);	
	fmpz_addmul_si(rank, product, combinations);
}

//SET_PART
//std::array<std::array< gen_rgf_table(int n, int k) {
//	
//}
/*
void rgf_rank(big rank, std::vector<uint8_t>& rgf, int k) {
	// This calculates the rank of a set-partition with exactly k-parts
	int n = rgf.size();
	int currentMax = 1;
	
	//Gen Tabl
	
	fmpz_zero(rank);
	
	// RGF always starts with 1, so we iterate from the second element
	for (int i = 1; i < n; i++) {
		int remLen = n - 1 - i;
		uint8_t digit = rgf[i];
		
		// We count the "branches" of the decision tree we skipped
		// Branches are ordered 1, 2, ..., currentMax, currentMax+1
		
		if (digit == currentMax+1) {
			// We picked the "new block" option.
			// This implies we skipped all 'currentMax' options to join existing blocks.
			// Each skipped option has weight table[remLen][currentMax]
			big countSkipped;
			fmpz_mul_ui(countSkipped, table[remLen][currentMax], currentMax);			
			fmpz_add(rank, rank, countSkipped);			
			currentMax++;
		}
		else {			
			// We picked an existing block (digit <= currentMax).
			// We skipped (digit - 1) options of joining smaller existing blocks.			
			fmpz_addmul_ui(rank, table[remLen][currentMax], digit-1);
			// currentMax does not change
		}
	}
	
}


std::vector<uint8_t> rgf_unrank(big rank, int n, int k) {
	//Converts a rank back into an RGF of length n with exactly k parts
	//Gen table
	std::vector<uint8_t> rgf(n);
	int currentMax = 1;
	
	for (int i = 1; i < n; i++) {
		int remLen = n - 1 - i;
		
		// Calculate the "weight" (number of possibilities) if we join an existing block
		big weightStay;
		weightStay = table[remLen][currentMax];
		
		// Total possibilities covered by joining ANY existing block (1..currentMax)
		big countStay;
		fmpz_mul_ui(countStay, weightStay, currentMax);
		
		if (fmpz_cmp(rank, countStay) < 0) {		
			// The rank falls within the "join existing block" range.
			// Determine exactly which block index (1..currentMax)
			big bigVal;
			fmpz_tdiv_q(bigVal, rank, weightStay
			int val = fmpz_get_ui(bigVal) + 1;
			rgf[i] = val;
			fmpz_mod(rank, rank, weightStay);
		}
		else {
			// The rank is beyond the "join existing" range, so we must start a new block.
			rgf[i] = currentMax + 1;
			fmpz_sub(rank, rank, countStay);			
			currentMax++;
		}
	}
	
	return rgf;
	
}
//END SET_PART
*/


void near_entropic_rank(big rank) {
			
	
	uint8_t valSeq[] = {1,7,7,1,14,7,0,11,2,13,13,12,11,2,0,7};
	if (DEBUG) {
		cout << "Ranking: ";
		printSeq(valSeq, SEQ_LEN);
	}
	
	//Pre-process file bits -> seq	
	int valToSym[MAX_SYM]; //Put val in get sym out
	for (int i = 0; i < MAX_SYM; i++) {
		valToSym[i] = INVALID;
	}
	std::vector<uint8_t> vals;
	uint8_t symSeq[SEQ_LEN];
		
	int symCount = INVALID;
	//Loop read in files as i16
	for (int i = 0; i < SEQ_LEN; i++){
		uint8_t val = valSeq[i];
		uint8_t sym;
		//Create map of seq -> vals
		if (valToSym[val] == INVALID) { //First time this symbol has been seen		
			symCount++;
			valToSym[val] = symCount;
			sym = symCount;			
			vals.push_back(val);
		}
		else sym = valToSym[val];
		
		symSeq[i] = sym;			
	}
	symCount++; //So that it reflects the total properly 
	
	//Permutation of symbols -> vals
	
	
	if (DEBUG) {
		cout << "Sym Count: " << symCount << endl;
		cout << "Sym Seq: ";
		printSeq(symSeq, SEQ_LEN);	
		
		cout << "Vals: ";
		printVector(vals);		
	}
	
	
	fmpz_zero(rank);	
	
	// 1. Add symbol sections
	for (int i = 1; i < symCount; i++) {				
		addSymbolSection(rank, i);
		//cout << "here" << endl;
	}
		
	if (DEBUG) {
		cout << "Rank after symbol section: " << endl;
		fmpz_print(rank);
		cout << endl;
	}
	
	//Calculate section sizes	
	uint64_t totalComb = comb(MAX_SYM, symCount);
	big combSectionSize;
	fmpz_fac_ui(combSectionSize, symCount);
	
	big stirSectionSize;
	fmpz_mul_si(stirSectionSize, combSectionSize, totalComb);
	
	
	//  2. Add Set Partition / Stirling2 rank 	
	//Note: the symSeq is just the 0-indexed version of the RGF
	big stirRank;
	//TODO!!!
	fmpz_addmul(rank, stirRank, stirSectionSize);
	
	if (DEBUG) {
		cout << "Stir Rank: " << endl;
		fmpz_print(stirRank);
		cout << endl;
	}
	
	//  3. Add the combination rank 
	uint64_t combRank = comb_rank(vals);	
	fmpz_addmul_si(rank, combSectionSize, combRank);
	if (DEBUG) {
		cout << "Comb Rank: " << combRank << endl;
		//std::vector<uint8_t> cv = {1, 2, 3, 8, 12, 13, 14, 15};
		//uint64_t cr = comb_rank(cv);
		//std::vector<uint8_t> ucr = comb_unrank(cr, MAX_SYM, 8);
		//cout << "Comb Rank: " << cr << endl;
		//cout << "Unrank: ";
		//printVector(ucr);
	}
	
	//  4. Add the Sym Perm Rank (Myrvold)	
	std::vector<uint8_t> symPerm(symCount);	
	for (int i = 0; i < MAX_SYM; i++) { 
		//Traverse vals in order
		if (valToSym[i] != INVALID) {
			symPerm.push_back(valToSym[i]);
		}
	}
	big myrvoldRank;
	myrvold_rank(myrvoldRank, symPerm);
	fmpz_add(rank, rank, myrvoldRank);
	
	if (DEBUG) {
		cout << "Sym Perm: ";
		printVector(symPerm);
		
		cout << "Myrvold Rank: ";
		fmpz_print (myrvoldRank);
		cout << endl;
	}
	
	if (DEBUG) {
		cout << "Final Rank: " << endl;
		fmpz_print(rank);
		cout << endl;
	}
}

int main() {
		
	big rank;
	near_entropic_rank( rank);	
	//fmpz_print(rank);

	return 0;
}