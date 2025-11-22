#include "permutations.h"

void myrvold_rank_recur(fmpz_t rank, int n, std::vector<uint8_t>& perm, std::vector<uint8_t>& invPerm) {
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
	myrvold_rank_recur(rank, n-1, perm, invPerm);
	fmpz_mul_ui(rank, rank, n);
	fmpz_add_ui(rank, rank, s);
	
}



void myrvold_rank(fmpz_t rank, std::vector<uint8_t>& perm) {
	fmpz_zero(rank);
	int permSize = perm.size();
	std::vector<uint8_t> invPerm(perm.size());
	for (int i = 0; i < permSize; i++) {
		uint8_t pos = perm[i];
		invPerm[pos] = i;
	}
	
	myrvold_rank_recur(rank, permSize, perm, invPerm);
}

void myrvold_unrank_recur(fmpz_t rank, int n, std::vector<uint8_t>& perm) {
	if (n < 1) return; //Anchor
	
	fmpz_t q;
	fmpz_tdiv_q_ui(q, rank, n);
	fmpz_t m;
	fmpz_mod_ui(m, rank, n);	
	int r = fmpz_get_ui(m);	
	
	//Swap
	uint8_t tmp1 = perm[n-1];
	uint8_t tmp2 = perm[r];
	perm[r] = tmp1;
	perm[n-1] = tmp2;
	
	myrvold_unrank_recur(q, n-1, perm);
	
}

void myrvold_unrank(std::vector<uint8_t>& perm, fmpz_t rank, int permSize) {	
	//Start with identity so that unranking can shuffle into place		
	
	for (int i = 0; i < permSize; i++) {
		perm[i] = i;
	}
		
	myrvold_unrank_recur(rank, permSize, perm); //Perm is mutated in the recursive function	
}