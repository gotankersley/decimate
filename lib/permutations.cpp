#include "permutations.h"
#include <iostream>

void myrvold_rank_recur(int n, std::vector<uint8_t>& perm, std::vector<uint8_t>& invPerm, fmpz_t rankOut) {
	if (n < 2) return; //Anchor
	
	uint8_t s = perm[n-1];
	//Swap
	uint8_t tmp1 = perm[invPerm[n-1]];
	uint8_t tmp2 = perm[n-1];
	perm[n-1] = tmp1;
	perm[invPerm[n-1]] = tmp2;
	
	tmp1 = invPerm[n-1];
	tmp2 = invPerm[s];
	invPerm[s] = tmp1;
	invPerm[n-1] = tmp2;
	
	//Recurse
	myrvold_rank_recur(n-1, perm, invPerm, rankOut);
	fmpz_mul_ui(rankOut, rankOut, n);
	fmpz_add_ui(rankOut, rankOut, s);
	
}



void myrvold_rank(std::vector<uint8_t> perm, fmpz_t rankOut) {
	fmpz_zero(rankOut);
	int permSize = perm.size();
	std::vector<uint8_t> invPerm(perm.size());
	for (int i = 0; i < permSize; i++) {
		uint8_t pos = perm[i];
		invPerm[pos] = i;
	}
	
	myrvold_rank_recur(permSize, perm, invPerm, rankOut);
}

void myrvold_unrank_recur(fmpz_t rank, int n, std::vector<uint8_t>& permOut) {
	if (n < 1) return; //Anchor
		
		
	fmpz_t m;
	fmpz_init(m);
	fmpz_mod_ui(m, rank, n);	
	int r = fmpz_get_ui(m);	
	fmpz_clear(m);
	
	fmpz_tdiv_q_ui(rank, rank, n);
	
	//Swap
	uint8_t tmp1 = permOut[n-1];
	uint8_t tmp2 = permOut[r];
	permOut[r] = tmp1;
	permOut[n-1] = tmp2;
	
	myrvold_unrank_recur(rank, n-1, permOut);
	
}

void myrvold_unrank(fmpz_t rank, std::vector<uint8_t>& permOut) {	
	//Start with identity so that unranking can shuffle into place		
	int permSize = permOut.size();
	for (int i = 0; i < permSize; i++) {
		permOut[i] = i;		
	}
		
	myrvold_unrank_recur(rank, permSize, permOut); //Perm is mutated in the recursive function	
}