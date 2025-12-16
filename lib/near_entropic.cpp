#include "near_entropic.h"
using std::cout, std::endl;


const bool DEBUG = false;
const int INVALID = -1;

//Combinatorial function to count the ways to rank:
// - Filling up a sequence of N length (e.g. N = 3)
// - Using an alphabet of A symbols available (e.g. A = {A,B,C} = 3 )
// - Using subset K of symbols (e.g. K = 2)
void addSymbolSection(fmpz_t rank, int seqLen, int maxSym, int k) {
	if (k < 0 || k > seqLen || k > maxSym) return;
	
	uint64_t combinations = comb(maxSym, k);
	
	fmpz_t factorialK;	
	fmpz_init(factorialK);
	fmpz_fac_ui(factorialK, k);	
	
	fmpz_t stir2;
	fmpz_init(stir2);
	arith_stirling_number_2(stir2, seqLen, k);

	fmpz_t product;
	fmpz_init(product);
	fmpz_mul(product, factorialK, stir2);	
	fmpz_addmul_ui(rank, product, combinations);
	
	fmpz_clear(factorialK);
	fmpz_clear(stir2);
	fmpz_clear(product);
}



void near_entropic_rank(std::vector<uint8_t>& valSeq, int maxSym, fmpz_t rankOut) {
			
		
	
	//Pre-process seq	
	int seqLen = valSeq.size();
		
	std::vector<int> valToSym(maxSym); //Put val in get sym out
	
	for (int i = 0; i < maxSym; i++) {
		valToSym[i] = INVALID;
	}

	std::vector<uint8_t> rgfSeq(seqLen);
			
	int symCount = INVALID;	 
	for (int i = 0; i < seqLen; i++){
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
	
	if (DEBUG) {
		cout << "Ranking: ";		
		printVector(valSeq);		
	}
	
	//Permutation of symbols -> vals
	std::vector<uint8_t> vals;
	std::vector<uint8_t> symPerm;	
	
	for (int i = 0; i < maxSym; i++) { 
		//Traverse vals in order
		if (valToSym[i] != INVALID) {
			vals.push_back(i);			
			symPerm.push_back(valToSym[i]);
		}		
	}
	
	if (DEBUG) {
		cout << "Sym Count: " << symCount << endl;
		cout << "RGF Seq: ";
		printVector(rgfSeq);
		
		cout << "Vals: ";
		printVector(vals);		
	}
	
	
	fmpz_zero(rankOut);	
	
	// 1. Add symbol sections
	for (int i = 1; i < symCount; i++) {				
		addSymbolSection(rankOut, seqLen, maxSym, i);		
	}
		
	if (DEBUG) {
		cout << "Rank after symbol section: " << endl;
		fmpz_print(rankOut);
		cout << endl;
	}
	
	//Calculate section sizes	
	uint64_t totalComb = comb(maxSym, symCount);
	fmpz_t combSectionSize;
	fmpz_init(combSectionSize);
	fmpz_fac_ui(combSectionSize, symCount);
	
	fmpz_t stirSectionSize;
	fmpz_init(stirSectionSize);
	fmpz_mul_ui(stirSectionSize, combSectionSize, totalComb);
	
	if (DEBUG) {
		cout << "Comb Section Size: " << endl;
		fmpz_print(combSectionSize);
		cout << endl;
		cout << "Stir Section Size: " << endl;
		fmpz_print(stirSectionSize);
		cout << endl;
	}
	
	//  2. Add Set Partition / Stirling2 rank 		
	fmpz_t stirRank;
	fmpz_init(stirRank);		
	rgf_rank(rgfSeq, symCount, stirRank);
	fmpz_addmul(rankOut, stirRank, stirSectionSize);
	fmpz_clear(stirSectionSize);
	
	if (DEBUG) {
		cout << "Stir Rank: " << endl;
		fmpz_print(stirRank);
		cout << endl;
	}
	fmpz_clear(stirRank);
	
	//  3. Add the combination rank 
	uint64_t combRank = comb_rank(vals);	
	fmpz_addmul_ui(rankOut, combSectionSize, combRank);
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
	if (DEBUG) {
		cout << "Sym Perm: ";
		printVector(symPerm);
	}
	myrvold_rank(symPerm, symRank);
	fmpz_add(rankOut, rankOut, symRank);
	
	if (DEBUG) {			
		cout << "Sym Rank: ";
		fmpz_print (symRank);
		cout << endl;
		
		cout << "Final Rank: " << endl;
		fmpz_print(rankOut);
		cout << endl;
	}
	fmpz_clear(symRank);			
}

void near_entropic_unrank(fmpz_t rank, int seqLen, int maxSym, std::vector<int>& countsOut, std::vector<uint8_t>& rgfOut) {
			
	
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
	for (int i = 1; i <= maxSym; i++) {
		addSymbolSection(count, seqLen, maxSym, i);		
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
	uint64_t totalComb = comb(maxSym, symCount);
	fmpz_t combSectionSize;
	fmpz_init(combSectionSize);
	fmpz_fac_ui(combSectionSize, symCount);
	
	fmpz_t stirSectionSize;
	fmpz_init(stirSectionSize);
	fmpz_mul_ui(stirSectionSize, combSectionSize, totalComb);
	

	if (DEBUG) {
		cout << "Comb Section Size: " << endl;
		fmpz_print(combSectionSize);
		cout << endl;
		cout << "Stir Section Size: " << endl;
		fmpz_print(stirSectionSize);
		cout << endl;
	}

	// 2. Get the values from the combination rank of symbols
	fmpz_t rankModStir;	
	fmpz_init(rankModStir);	
	fmpz_mod(rankModStir, rank, stirSectionSize);			
	fmpz_tdiv_q(rankModStir, rankModStir, combSectionSize);	
	uint64_t combRank = fmpz_get_ui(rankModStir);	
	fmpz_clear(rankModStir);	
	std::vector<uint8_t> combVals(symCount);	
	comb_unrank(combRank, maxSym, symCount, combVals);
	
	
	if (DEBUG) {
		cout << "Comb Rank: " << combRank << endl;
		cout << "Comb Vals: ";
		printVector(combVals);
	}
	
	//  3. Get the Sym Perm from Myrvold Rank	
	fmpz_t symRank;
	fmpz_init(symRank);
	fmpz_mod(symRank, rank, combSectionSize);	
	std::vector<uint8_t> symPerm(symCount);
	if (DEBUG) {	
		cout << "Sym Rank: ";
		fmpz_print (symRank);
		cout << endl;
	}
	myrvold_unrank(symRank, symPerm);
	fmpz_clear(symRank);
	fmpz_clear(combSectionSize);
	std::vector<uint8_t> invPerm(symCount);
	for (int i = 0; i < symCount; i++) {
		invPerm[symPerm[i]] = i;
	}
	
	if (DEBUG) {
		cout << "Sym Perm: ";
		printVector(symPerm);		
	}

	// 4. Get the Set Partition from Stirling2 rank		
	// Note: We're also applying inverse perm to recreate seq	
	fmpz_t stirRank;
	fmpz_init(stirRank);
	fmpz_tdiv_q(stirRank, rank, stirSectionSize);
	fmpz_clear(stirSectionSize);
	if (DEBUG) {
		cout << "Stir Rank: " << endl;
		fmpz_print(stirRank);
		cout << endl;					
	}	
	rgf_unrank_opt(stirRank, seqLen, symCount, combVals, invPerm, countsOut, rgfOut);		
	fmpz_clear(stirRank);
	if (DEBUG) {
		cout << "RGF Seq: ";
		printVector(rgfOut);	
	}
		
		
	
	if (DEBUG) {
		//cout << "Entropy of Unranked Seq: "<< std::fixed << std::setprecision(2) << measureEntropy(counts, seqLen) << endl;
		cout << "Final Seq: ";
		printVector(rgfOut);
	}
}



