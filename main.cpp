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
		cout << +val << ","; 
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
	fmpz_init(factorialK);
	fmpz_fac_ui(factorialK, k);	
	
	fmpz_t stir2;
	fmpz_init(stir2);
	arith_stirling_number_2(stir2, SEQ_LEN, k);

	fmpz_t product;
	fmpz_init(product);
	fmpz_mul(product, factorialK, stir2);	
	fmpz_addmul_si(rank, product, combinations);
	
	fmpz_clear(factorialK);
	fmpz_clear(stir2);
	fmpz_clear(product);
}



void near_entropic_rank(fmpz_t rank, std::vector<uint8_t>& valSeq) {
			
	
	//uint8_t valSeq[] = {1,7,7,1,14,7,0,11,2,13,13,12,11,2,0,7};
	if (DEBUG) {
		cout << "Ranking: ";		
		printVector(valSeq);
	}
	
	//Pre-process file bits -> seq	
	int valToSym[MAX_SYM]; //Put val in get sym out
	for (int i = 0; i < MAX_SYM; i++) {
		valToSym[i] = INVALID;
	}

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
		}
		else sym = valToSym[val];
		
		rgfSeq[i] = sym + 1; //RGF is 1-indexed			
	}
	symCount++; //So that it reflects the total properly 
	
	//Permutation of symbols -> vals
	std::vector<uint8_t> vals;
	std::vector<uint8_t> symPerm(symCount);	
	for (int i = 0; i < MAX_SYM; i++) { 
		//Traverse vals in order
		if (valToSym[i] != INVALID) {
			vals.push_back(i);
			symPerm[i] = valToSym[i];
		}		
	}
	
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
	fmpz_init(combSectionSize);
	fmpz_fac_ui(combSectionSize, symCount);
	
	fmpz_t stirSectionSize;
	fmpz_init(stirSectionSize);
	fmpz_mul_si(stirSectionSize, combSectionSize, totalComb);
	
	//  2. Add Set Partition / Stirling2 rank 		
	fmpz_t stirRank;
	fmpz_init(stirRank);
	rgf_rank(rgfSeq, symCount, stirRank);
	fmpz_addmul(rank, stirRank, stirSectionSize);
	fmpz_clear(stirSectionSize);
	
	if (DEBUG) {
		cout << "Stir Rank: " << endl;
		fmpz_print(stirRank);
		cout << endl;
	}
	fmpz_clear(stirRank);
	
	//  3. Add the combination rank 
	uint64_t combRank = comb_rank(vals);	
	fmpz_addmul_si(rank, combSectionSize, combRank);
	if (DEBUG) {
		cout << "Comb Rank: " << combRank << endl;		
		cout << "Comb Vals: ";
		printVector(vals);
		cout << "Comb Section Size: ";
		fmpz_print(combSectionSize);
		cout << endl;
	}
	fmpz_clear(combSectionSize);
	
	//  4. Add the Sym Perm Rank (Myrvold)	
	fmpz_t symRank;
	fmpz_init(symRank);
	myrvold_rank(symPerm, symRank);
	fmpz_add(rank, rank, symRank);
	
	if (DEBUG) {
		cout << "Sym Perm: ";
		printVector(symPerm);
		
		cout << "Sym Rank: ";
		fmpz_print (symRank);
		cout << endl;
		
		cout << "Final Rank: " << endl;
		fmpz_print(rank);
		cout << endl;
	}
	fmpz_clear(symRank);			
}

void near_entropic_unrank(fmpz_t rank, std::vector<uint8_t>& rgfSeq) {
			
	
	//  1. Get symbol section
	if (DEBUG) {
		cout << "Unranking" << endl;
		cout << "Rank: " << endl;
		fmpz_print(rank);
		cout << endl;
	}
	fmpz_t count;
	fmpz_t prevCount;
	
	fmpz_init(count);
	fmpz_init(prevCount);
	
	fmpz_zero(count);
	fmpz_zero(prevCount);
	int symCount = INVALID;
	for (int i = 1; i <= MAX_SYM; i++) {
		addSymbolSection(count, i);		
		if (fmpz_cmp(count, rank) > 0) {		
			symCount = i;
			fmpz_sub(rank, rank, prevCount);
			break;
		}
		else fmpz_set(prevCount, count);
	}
	
	if (DEBUG) {
		cout << "Sym Count: " << symCount << endl;		
	}
	fmpz_clear(count);
	fmpz_clear(prevCount);
	
	
	//Calculate section sizes	
	uint64_t totalComb = comb(MAX_SYM, symCount);
	fmpz_t combSectionSize;
	fmpz_init(combSectionSize);
	fmpz_fac_ui(combSectionSize, symCount);
	
	fmpz_t stirSectionSize;
	fmpz_init(stirSectionSize);
	fmpz_mul_si(stirSectionSize, combSectionSize, totalComb);
	

	// 2. Get the Set Partition from Stirling2 rank		
	fmpz_t stirRank;
	fmpz_init(stirRank);
	fmpz_tdiv_q(stirRank, rank, stirSectionSize);
	
	rgf_unrank(stirRank, SEQ_LEN, symCount, rgfSeq);
	fmpz_addmul(rank, stirRank, stirSectionSize);
	
	if (DEBUG) {
		cout << "Stir Rank: " << endl;
		fmpz_print(stirRank);
		cout << endl;
	}
	fmpz_clear(stirRank);
	
	// 3. Get the values from the combination rank of symbols
	fmpz_t rankModStir;
	fmpz_init(rankModStir);
	fmpz_mod(rankModStir, rank, stirSectionSize);	
	fmpz_clear(stirSectionSize);
	fmpz_tdiv_q(rankModStir, rankModStir, combSectionSize);
	uint64_t combRank = fmpz_get_ui(rankModStir);
	fmpz_clear(rankModStir);
	std::vector<uint8_t> combVals;
	comb_unrank(combRank, MAX_SYM, symCount, combVals);
	
	if (DEBUG) {
		cout << "Comb Rank: " << combRank << endl;
	}
	
	//  4. Add the Sym Perm Rank (Myrvold)	
	fmpz_t symRank;
	fmpz_init(symRank);
	fmpz_mod(symRank, rank, combSectionSize);	
	std::vector<uint8_t> symPerm(symCount);
	myrvold_unrank(symRank, symPerm);
	fmpz_clear(combSectionSize);
	if (DEBUG) {
		cout << "Sym Perm: ";
		printVector(symPerm);
		
		cout << "Myrvold Rank: ";
		fmpz_print (symRank);
		cout << endl;
	}
	fmpz_clear(symRank);
	
	// Apply to recreate seq
	std::vector<uint8_t> symToVal(symCount);
	for (int i = 0; i < (int)combVals.size(); i++) {
		uint8_t val = combVals[i];
		symToVal[i] = val;
	}
	//TODO - this could maybe be moved to the Stirling section to avoid another loop over SEQ_LEN
	for (int i = 0; i < SEQ_LEN; i++) {
		uint8_t val = symToVal[symPerm[(rgfSeq[i]-1)]];
		rgfSeq[i] = val;
	}
	
	if (DEBUG) {
		cout << "Final Seq: ";
		printVector(rgfSeq);
	}
}

int main() {
		
	std::vector<uint8_t> valSeq = {1,7,7,1,14,7,0,11,2,13,13,12,11,2,0,7};
	fmpz_t rank;
	fmpz_init(rank);
	near_entropic_rank(rank, valSeq);			
	//std::vector<uint8_t> rgfSeq(SEQ_LEN);
	//near_entropic_unrank(rank, rgfSeq);	
	fmpz_clear(rank);
	return 0;
}