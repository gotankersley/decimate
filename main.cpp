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


void printSeq(const sym_t arr[], int size) {    
	cout << "[";
    for (int i = 0; i < size; ++i) {
        cout << +arr[i] << ",";
    }
    cout << "]" << endl;
}

void printVector(std::vector<sym_t> vals) {    
	for (int val : vals) {
		cout << val << ","; 
	}
	cout << endl;
}

//COMBO LIB
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

uint64_t comb_rank(std::vector<sym_t> vals) {
	int k = vals.size();
	uint64_t rank = 0;
	
	for (int i = 0; i < k; i++) {
		rank += comb(vals[k-1-i], k-i);
	}
	return rank;
}

std::vector<sym_t> comb_unrank(uint64_t rank, int n, int k) {
	std::vector<sym_t> vals(k);
	
	for (int i = 0; i < k; i++) {
		while (comb(n, k-i) > rank) {
			n--;
		}
		vals[k-i-1] = n;
		rank -= comb(n, k - i);
	}
	return vals;
}

//END COMBO LIB

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




void near_entropic_rank(big rank) {
	std::vector<sym_t> cv = {1, 2, 3, 8, 12, 13, 14, 15};
	uint64_t cr = comb_rank(cv);
	std::vector<sym_t> ucr = comb_unrank(cr, MAX_SYM, 8);
	cout << "Comb Rank: " << cr << endl;
	cout << "Unrank: ";
	printVector(ucr);
	
	sym_t valSeq[] = {1,7,7,1,14,7,0,11,2,13,13,12,11,2,0,7};
	if (DEBUG) {
		cout << "Ranking: ";
		printSeq(valSeq, SEQ_LEN);
	}
	
	//Pre-process file bits -> seq	
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
	}
	
	//  3. Add the combination rank 
	
	
	//  4. Add the Sym/Set-Part Perm Rank (Myrvold)	
	
	if (DEBUG) {
		cout << "Final Rank: " << endl;
		fmpz_print(rank);
	}
}

int main() {
		
	big rank;
	near_entropic_rank( rank);	
	//fmpz_print(rank);

	return 0;
}