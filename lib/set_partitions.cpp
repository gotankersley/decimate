#include "set_partitions.h"

void gen_rgf_table(fmpz_mat_t table, int n, int k) {
	//Generates a table where table[remLen][currentMax] stores the number of ways 
	//to complete a partition of 'remLen' remaining elements, given the 
	//'currentMax' block index, such that the final number of blocks is exactly k
	
	// Table dimensions: (n+1) rows for remaining length, (k+2) cols for current max
	fmpz_mat_init(table, n+1, k+2);
	
	// Base Case: If 0 elements remain, we represent a valid partition ONLY if
	// we have already reached exactly k blocks	
	fmpz_one(fmpz_mat_entry(table, 0, k));
	
	for (int length = 1; length < (n+1); length++) {
		for (int m = 1; m < (k+1); m++) {
			// Recurrence:
			// 1. Join an existing block (1..m). There are 'm' choices.
			//	 We consume 1 length, max stays 'm'.
			// 2. Create a new block (m+1). There is 1 choice.
			//	 We consume 1 length, max becomes 'm+1'.
			// Note: If m+1 > k, table[length-1][m+1] will be 0, enforcing the limit.
			
			//table[length][m] = (m * table[length-1][m]) + table[length-1][m+1]	
			fmpz_t product;
			fmpz_mul_ui(product, fmpz_mat_entry(table, length-1, m), m);
			fmpz_add(
				fmpz_mat_entry(table, length, m),
				product,
				fmpz_mat_entry(table, length-1, m+1)
			);
		}
	}
	
	fmpz_mat_print_pretty(table);
	//Don't forget to call fmpz_mat_clear(table);
}

void rgf_rank(fmpz_t rank, std::vector<uint8_t>& rgf, int k) {
	// This calculates the rank of a set-partition with exactly k-parts
	int n = rgf.size();
	int currentMax = 1;
	
	fmpz_mat_t table;
	gen_rgf_table(table, n, k);
	
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
			fmpz_t countSkipped;
			fmpz_mul_ui(countSkipped, fmpz_mat_entry(table, remLen, currentMax), currentMax);			
			fmpz_add(rank, rank, countSkipped);			
			currentMax++;
		}
		else {			
			// We picked an existing block (digit <= currentMax).
			// We skipped (digit - 1) options of joining smaller existing blocks.			
			fmpz_addmul_ui(rank, fmpz_mat_entry(table, remLen, currentMax), digit-1);
			// currentMax does not change
		}
	}
	fmpz_mat_clear(table);
}


void rgf_unrank(std::vector<uint8_t>& rgf, fmpz_t rank, int n, int k) {
	//Converts a rank back into an RGF of length n with exactly k parts
	//std::vector<uint8_t> rgf(n);
	fmpz_mat_t table;
	gen_rgf_table(table, n, k);
	
	int currentMax = 1;
	
	for (int i = 1; i < n; i++) {
		int remLen = n - 1 - i;
		
		// Calculate the "weight" (number of possibilities) if we join an existing block
		fmpz_t weightStay;
		fmpz_set(weightStay, fmpz_mat_entry(table, remLen, currentMax));
		
		// Total possibilities covered by joining ANY existing block (1..currentMax)
		fmpz_t countStay;
		fmpz_mul_ui(countStay, weightStay, currentMax);
		
		if (fmpz_cmp(rank, countStay) < 0) {		
			// The rank falls within the "join existing block" range.
			// Determine exactly which block index (1..currentMax)
			fmpz_t fmpz_tVal;
			fmpz_tdiv_q(fmpz_tVal, rank, weightStay);
			int val = fmpz_get_ui(fmpz_tVal) + 1;
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
	fmpz_mat_clear(table);		
}