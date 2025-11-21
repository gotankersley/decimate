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
		for (int val : vals) {
			cout << val << ","; 
		}
		cout << endl;
	}
	
	
	fmpz_zero(rank);	
	
	// 1. Add symbol sections
	for (int i = 1; i < symCount; i++) {				
		addSymbolSection(rank, i);
		//cout << "here" << endl;
	}
		
	//if (DEBUG) count << "rank after symbol section: <<  rank << endl;
	
	//arith_stirling_number_2(rank, n, k);	
}

int main() {
		
	big rank;
	near_entropic_rank( rank);	
	//fmpz_print(rank);

	return 0;
}