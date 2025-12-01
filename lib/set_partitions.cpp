#include "set_partitions.h"
#include <iostream>

void gen_rgf_table(int n, int k, fmpz_mat_t tableOut) {
	//Generates a table where table[remLen][currentMax] stores the number of ways 
	//to complete a partition of 'remLen' remaining elements, given the 
	//'currentMax' block index, such that the final number of blocks is exactly k
	
	// Table dimensions: (n+1) rows for remaining length, (k+2) cols for current max
	fmpz_mat_init(tableOut, n+1, k+2);
	
	// Base Case: If 0 elements remain, we represent a valid partition ONLY if
	// we have already reached exactly k blocks	
	fmpz_one(fmpz_mat_entry(tableOut, 0, k));
	
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
			fmpz_mul_ui(product, fmpz_mat_entry(tableOut, length-1, m), m);
			fmpz_add(
				fmpz_mat_entry(tableOut, length, m),
				product,
				fmpz_mat_entry(tableOut, length-1, m+1)
			);
		}
	}
	
	//fmpz_mat_print_pretty(tableOut);
	//Don't forget to call fmpz_mat_clear(table);
}

void gen_rgf_row_old(int n, int k, uint8_t& CUR, fmpz_mat_t rowOut) { //Deprecated
	//Computes only the N-th row of partial sums of the Stirling table
    
	fmpz_mat_init(rowOut, 2, k+2); //Two rows to swap to optimize memory allocation
	fmpz_one(fmpz_mat_entry(rowOut, 0, k)); // Base case: length 0 is valid only if max is already k
	
	uint8_t NEXT = !CUR;
	
	for (int length = 1; length < (n+1); length++) {		
		if (length % 1000 == 0) std::cout << length << std::endl;
		for (int m = 1; m < (k+1); m++) {
			// Standard Recurrence: S(L, m) = m*S(L-1, m) + S(L-1, m+1)
			fmpz_t product;
			fmpz_mul_ui(product, fmpz_mat_entry(rowOut, CUR, m), m);
			fmpz_add(
				fmpz_mat_entry(rowOut, NEXT, m),
				product,
				fmpz_mat_entry(rowOut, CUR, m+1)
			);
		}
		//Efficient bitwise swap
		CUR ^= NEXT;
		NEXT ^= CUR;
		CUR ^= NEXT;
	}	
}

void gen_rgf_row(int n, int k, fmpz_mat_t rowOut) {
	//Note: This calculates an arbitrary starting row, which is signicantly quicker
	//than calculating from the beginning
	fmpz_mat_init(rowOut, 2, k+2); //Two rows to swap to optimize memory allocation	
			
	for (int m = 1; m < (k+1); m++) {
		int level = m-1;
		gen_rgf_cell(n, k, level, rowOut);	
	}
}

void gen_rgf_cell(int n, int k, int level, fmpz_mat_t rowOut) {
	//Calculate arbitrary cell of RGF table
	//Note: This can be used to calculate a starting row significantly
	//quicker than calculating the whole table.  However, if the entire table
	//needs to be used, (as in the case of ranking), the incremental approach probably is better.		
	int col = k-level;
	uint8_t CUR = 0;
	fmpz_zero(fmpz_mat_entry(rowOut, CUR, col));
	if (level >= k) return;
		
	
	fmpz_t pascalWeight; //Treating Pascal Triangle like a lookup table of weights
	fmpz_t power;
	
	fmpz_init (pascalWeight);
	fmpz_init (power);
	
	for (int i = 0; i < level+1; i++) {
		
		fmpz_bin_uiui(pascalWeight, level, i);
		
		if (i % 2 != 0) {
			fmpz_mul_si(pascalWeight, pascalWeight, -1);
		}
		
		fmpz_ui_pow_ui(power, k-i, n);
		fmpz_addmul(fmpz_mat_entry(rowOut, CUR, col), pascalWeight, power);
		
	}
	fmpz_clear(pascalWeight);
	fmpz_clear(power);
	
	fmpz_t factLevel;
	fmpz_init(factLevel);
	fmpz_fac_ui(factLevel, level);	
	fmpz_tdiv_q(
		fmpz_mat_entry(rowOut, CUR, col), 
		fmpz_mat_entry(rowOut, CUR, col), 
		factLevel
	);
	fmpz_clear(factLevel);
}

/*
void serialize_mat(const char* filename, fmpz_mat_t mat) {
	FILE* file = fopen(filename, "w");
    if (!file) {
        perror("Could not open file for writing");
        return;
    }
	
	long rows = mat->r;
	long cols = mat->c;	
	
	//Write dimensions
	fwrite(&rows, sizeof(long), 1, file);
	fwrite(&cols, sizeof(long), 1, file);
	
	//Write each element in raw binary format
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {			
			fmpz_out_raw(file, fmpz_mat_entry( mat, i, j));			
		}
	}
	fclose(file);	
}

void deserialize_mat(const char* filename, fmpz_mat_t mat) {
	FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Could not open file for reading");
        return;
    }
	
	long rows = mat->r;
	long cols = mat->c;	

	//Read dimensions
	fread(&rows, sizeof(long), 1, file);
	fread(&cols, sizeof(long), 1, file);
	
	//Init matrix
	fmpz_mat_init(mat, rows, cols);
	
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {		
			fmpz_inp_raw(fmpz_mat_entry( mat, i, j), file);			
		}
	}
	
	fclose(file);
	//fmpz_mat_print_pretty(mat);
}
*/
void rgf_rank(std::vector<uint8_t>& rgf, int k, fmpz_t rankOut) {
		
	// Get the starting row:
    // At index i=1, the remaining length is (n-2).
    // We compute the row for length = n-2.
	fmpz_mat_t row;
	gen_rgf_row(rgf.size()-2, k, row);	
	rgf_rank_row(rgf, k, row, rankOut);
	fmpz_mat_clear(row);
}


void rgf_rank_table(std::vector<uint8_t>& rgf, int k, fmpz_mat_t table, fmpz_t rankOut) {
	// This calculates the rank of a set-partition with exactly k-parts
	int n = rgf.size();
	int currentMax = 1;
	
	fmpz_zero(rankOut);
	
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
			fmpz_init(countSkipped);
			fmpz_mul_ui(countSkipped, fmpz_mat_entry(table, remLen, currentMax), currentMax);			
			fmpz_add(rankOut, rankOut, countSkipped);
			fmpz_clear(countSkipped);		
			currentMax++;
		}
		else {			
			// We picked an existing block (digit <= currentMax).
			// We skipped (digit - 1) options of joining smaller existing blocks.			
			fmpz_addmul_ui(rankOut, fmpz_mat_entry(table, remLen, currentMax), digit-1);
			// currentMax does not change
		}
	}	
}


void rgf_rank_row(std::vector<uint8_t>& rgf, int k, fmpz_mat_t row, fmpz_t rankOut) {
	//Ranks an RGF forward (index 1 -> n) using only O(k) space
    //by mathematically inverting the Stirling recurrence at each step.
	//NOTE: row is actually two rows - current and previous and we swap 
	//between them for more efficient memory usage
	int n = rgf.size();
	fmpz_zero(rankOut);
	if (n <= 1) return;
	
	
	uint8_t CUR = 0;
	uint8_t PREV = 1;
    
	int currentMax = 1;
	
	for (int i = 1; i < n; i++) {
		uint8_t digit = rgf[i];
		
		// Standard Forward Ranking Logic
        if (digit == currentMax + 1) {
            // We picked a new block. 
            // We skipped 'currentMax' branches that would have kept the max same.            
			fmpz_addmul_ui(rankOut, fmpz_mat_entry(row, CUR, currentMax), currentMax);
            currentMax++;
		}
		else {
			// We picked an existing block.
            // We skipped (digit - 1) branches of smaller existing blocks.
			fmpz_addmul_ui(rankOut, fmpz_mat_entry(row, CUR, currentMax), (digit - 1));
		}
		
		// Downgrade the table row for the next step (Length L -> L-1)
        // We only need to do this if there are steps remaining
		if (i < n - 1) {
			// Inverted Recurrence: S(L-1, m) = (S(L, m) - S(L-1, m+1)) / m
            // We must iterate backwards (k -> 1) because we need prev_row[m+1]
            for (int m = k; m >= currentMax; m--) {

                // Always divides evenly
				fmpz_sub(
					fmpz_mat_entry(row, PREV, m),				
					fmpz_mat_entry(row, CUR, m),				
					fmpz_mat_entry(row, PREV, m+1)
				);
				fmpz_tdiv_q_ui(
					fmpz_mat_entry(row, PREV, m), 
					fmpz_mat_entry(row, PREV, m),
					m
				);  
				//Note: The ranking can potentially be done in reverse,
				//which allows the minor optimization of multiplication, instead of division.
				//However, this is ONLY for ranking, as the unranking can not be done in reverse
			}
			
			//Swap
			CUR ^= PREV;
			PREV ^= CUR;
			CUR ^= PREV;			
		}
	}	
}

void rgf_unrank(fmpz_t rank, int n, int k, std::vector<uint8_t>& rgfOut) {		
	fmpz_mat_t row;	
	gen_rgf_row(n-1, k, row);
	rgf_unrank_row(rank, n, k, row, rgfOut);
	//std::cout << "Finished unrank gen" << std::endl;
	fmpz_mat_clear(row);
}
void rgf_unrank_opt(fmpz_t rank, int n, int k, std::vector<uint8_t>& combVals, std::vector<uint8_t>& invPerm, std::vector<int>& countsOut, std::vector<uint8_t>& rgfOut) {
	fmpz_mat_t row;	
	gen_rgf_row(n-1, k, row);
	rgf_unrank_row_opt(rank, n, k, row, combVals, invPerm, countsOut, rgfOut);	
	fmpz_mat_clear(row);
}
void rgf_unrank_table(fmpz_t rank, int n, int k, fmpz_mat_t table, std::vector<uint8_t>& rgfOut) {
	//Converts a rank back into an RGF of length n with exactly k parts	

	rgfOut[0] = 1;
	int currentMax = 1;
	
	fmpz_t countStay;
	fmpz_t tVal;
	fmpz_init(countStay);
	fmpz_init(tVal);
	for (int i = 1; i < n; i++) {
		int remLen = n - 1 - i;
		
		// Calculate the "weight" (number of possibilities) if we join an existing block
		fmpz_t weightStay;		
		fmpz_set(weightStay, fmpz_mat_entry(table, remLen, currentMax));
		
		// Total possibilities covered by joining ANY existing block (1..currentMax)		
		fmpz_mul_ui(countStay, weightStay, currentMax);
		
		if (fmpz_cmp(rank, countStay) < 0) {		
			// The rank falls within the "join existing block" range.
			// Determine exactly which block index (1..currentMax)
			
			fmpz_tdiv_q(tVal, rank, weightStay);
			int val = fmpz_get_ui(tVal) + 1;
			rgfOut[i] = val;
			fmpz_clear(tVal);
			fmpz_mod(rank, rank, weightStay);
		}
		else {
			// The rank is beyond the "join existing" range, so we must start a new block.
			rgfOut[i] = currentMax + 1;
			fmpz_sub(rank, rank, countStay);			
			currentMax++;
		}		
	}
	fmpz_clear(tVal);
	fmpz_clear(countStay);
	
}

void rgf_unrank_row(fmpz_t rank, int n, int k, fmpz_mat_t row, std::vector<uint8_t>& rgfOut) {
	//Unranks an RGF using O(K) space by inverting the Stirling recurrence on the fly.

	rgfOut[0] = 1;
	int currentMax = 1;
	uint8_t CUR = 0;
	uint8_t PREV = 1;
	
	fmpz_t countStay;
	fmpz_t tVal;
	fmpz_init(countStay);
	fmpz_init(tVal);
	for (int i = 1; i < n; i++) {
		if (i % 1000 == 0) std::cout << "Unrank: " << i << std::endl;
		// Inverted Recurrence: S(L-1, m) = (S(L, m) - S(L-1, m+1)) / m
		// We must iterate backwards (k -> 1) because we need prev_row[m+1]
		for (int m = k; m >= currentMax; m--) {

			// Always divides evenly
			fmpz_sub(
				fmpz_mat_entry(row, PREV, m),				
				fmpz_mat_entry(row, CUR, m),				
				fmpz_mat_entry(row, PREV, m+1)
			);
			fmpz_tdiv_q_ui(
				fmpz_mat_entry(row, PREV, m), 
				fmpz_mat_entry(row, PREV, m),
				m
			);                
		}
		
		//Swap
		CUR ^= PREV;
		PREV ^= CUR;
		CUR ^= PREV;

		// Now current row corresponds to the correct 'rem_len' for this part
		
		// Calculate the "weight" (number of possibilities) if we join an existing block
		fmpz_t weightStay;
		fmpz_set(weightStay, fmpz_mat_entry(row, CUR, currentMax));
		
		fmpz_mul_ui(countStay, weightStay, currentMax);
		
		if (fmpz_cmp(rank, countStay) < 0) {	
			// Stay with existing block
			fmpz_tdiv_q(tVal, rank, weightStay);
			int val = fmpz_get_ui(tVal) + 1;
			rgfOut[i] = val;
			fmpz_clear(tVal);
			fmpz_mod(rank, rank, weightStay);
		}
		else {
			// Create new block
			rgfOut[i] = currentMax + 1;
			fmpz_sub(rank, rank, countStay);			
			currentMax++;
		}		
		
	}
	fmpz_clear(tVal);
	fmpz_clear(countStay);	
}

void rgf_unrank_row_opt(fmpz_t rank, int n, int k, fmpz_mat_t row, std::vector<uint8_t>& combVals, std::vector<uint8_t>& invPerm, std::vector<int>& countsOut, std::vector<uint8_t>& rgfOut) {
	//Unranks an RGF using O(K) space by inverting the Stirling recurrence on the fly.

	rgfOut[0] = 1;
	int currentMax = 1;
	uint8_t CUR = 0;
	uint8_t PREV = 1;
	
	fmpz_t countStay;
	fmpz_t tVal;
	fmpz_init(countStay);
	fmpz_init(tVal);
	for (int i = 1; i < n; i++) {
		if (i % 1000 == 0) std::cout << "Unrank: " << i << std::endl;
		// Inverted Recurrence: S(L-1, m) = (S(L, m) - S(L-1, m+1)) / m
		// We must iterate backwards (k -> 1) because we need prev_row[m+1]
		for (int m = k; m >= currentMax; m--) {

			// Always divides evenly
			fmpz_sub(
				fmpz_mat_entry(row, PREV, m),				
				fmpz_mat_entry(row, CUR, m),				
				fmpz_mat_entry(row, PREV, m+1)
			);
			fmpz_tdiv_q_ui(
				fmpz_mat_entry(row, PREV, m), 
				fmpz_mat_entry(row, PREV, m),
				m
			);                
		}
		
		//Swap
		CUR ^= PREV;
		PREV ^= CUR;
		CUR ^= PREV;

		// Now current row corresponds to the correct 'rem_len' for this part
		
		// Calculate the "weight" (number of possibilities) if we join an existing block
		fmpz_t weightStay;
		fmpz_set(weightStay, fmpz_mat_entry(row, CUR, currentMax));
		
		fmpz_mul_ui(countStay, weightStay, currentMax);
		
		if (fmpz_cmp(rank, countStay) < 0) {	
			// Stay with existing block
			fmpz_tdiv_q(tVal, rank, weightStay);
			int val = fmpz_get_ui(tVal) + 1;
			rgfOut[i] = val;
			fmpz_clear(tVal);
			fmpz_mod(rank, rank, weightStay);
		}
		else {
			// Create new block
			rgfOut[i] = currentMax + 1;
			fmpz_sub(rank, rank, countStay);			
			currentMax++;
		}		
		
		//Optimization - apply values here so we don't have to do another loop over the sequence
		uint8_t sym = rgfOut[i]-1;
		uint8_t combVal = combVals[invPerm[sym]];
		rgfOut[i] = combVal;
		countsOut[sym]++;
	}
	fmpz_clear(tVal);
	fmpz_clear(countStay);	
	
	//First value special case
	uint8_t sym = rgfOut[0]-1;
	uint8_t val = combVals[invPerm[sym]];
	rgfOut[0] = val;
	countsOut[sym]++;
}