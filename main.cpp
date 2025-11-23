#include <iostream>
#include <cstdint>
#include <vector>


#include "flint/fmpz.h"
#include "flint/arith.h"
#include "lib/combinations.h"
#include "lib/permutations.h"
#include "lib/set_partitions.h"
using std::cout, std::endl;



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
void addSymbolSection(fmpz_t rank, int k) {
	if (k < 0 || k > SEQ_LEN || k > MAX_SYM) return;
	
	uint64_t combinations = comb(MAX_SYM, k);
	
	fmpz_t factorialK;
	fmpz_fac_ui(factorialK, k);
	
	fmpz_t stir2;
	arith_stirling_number_2(stir2, SEQ_LEN, k);

	fmpz_t product;
	fmpz_mul(product, factorialK, stir2);	
	fmpz_addmul_si(rank, product, combinations);
}



void near_entropic_rank(fmpz_t rank) {
			
	
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
	std::vector<uint8_t> rgfSeq(SEQ_LEN);
		
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
		
		rgfSeq[i] = sym + 1; //RGF is 1-indexed			
	}
	symCount++; //So that it reflects the total properly 
	
	//Permutation of symbols -> vals
	
	
	if (DEBUG) {
		cout << "Sym Count: " << symCount << endl;
		cout << "Sym Seq: ";
		printVector(rgfSeq);
		
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
	fmpz_t combSectionSize;
	fmpz_fac_ui(combSectionSize, symCount);
	
	fmpz_t stirSectionSize;
	fmpz_mul_si(stirSectionSize, combSectionSize, totalComb);
	
	
	//  2. Add Set Partition / Stirling2 rank 		
	fmpz_t stirRank;
	rgf_rank(stirRank, rgfSeq, symCount);
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
		//std::vector<uint8_t> vals(k);
		//std::vector<uint8_t> ucr = comb_unrank(vals, cr, MAX_SYM, 8);
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
	fmpz_t myrvoldRank;
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
		
	fmpz_t rank;
	near_entropic_rank( rank);	
	//fmpz_print(rank);

	return 0;
}