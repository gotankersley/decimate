#include "set_partitions.h"

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

void rgf_rank(std::vector<uint8_t>& rgf, int k, fmpz_t rankOut) {
	fmpz_mat_t table;
	//gen_rgf_table(rgf.size(), k, table);
	std::string filename = std::to_string(k) + ".mat";
	deserialize_mat(filename.c_str() , table);	

	rgf_rank_table(rgf, k, table, rankOut);
	fmpz_mat_clear(table);
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
			fmpz_mul_ui(countSkipped, fmpz_mat_entry(table, remLen, currentMax), currentMax);			
			fmpz_add(rankOut, rankOut, countSkipped);			
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


void rgf_unrank(fmpz_t rank, int n, int k, std::vector<uint8_t>& rgfOut) {
	fmpz_mat_t table;
	//gen_rgf_table(n, k, table);
	std::string filename = std::to_string(k) + ".mat";
	deserialize_mat(filename.c_str() , table);	
	
	rgf_unrank_table(rank, n, k, table, rgfOut);
	fmpz_mat_clear(table);
}
void rgf_unrank_table(fmpz_t rank, int n, int k, fmpz_mat_t table, std::vector<uint8_t>& rgfOut) {
	//Converts a rank back into an RGF of length n with exactly k parts	
	//fmpz_mat_t table;
	//gen_rgf_table( n, k, table);
	rgfOut[0] = 1;
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
			rgfOut[i] = val;
			fmpz_mod(rank, rank, weightStay);
		}
		else {
			// The rank is beyond the "join existing" range, so we must start a new block.
			rgfOut[i] = currentMax + 1;
			fmpz_sub(rank, rank, countStay);			
			currentMax++;
		}
	}
	//fmpz_mat_clear(table);		
}