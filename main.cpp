#include <iostream>
#include <cstdint>
#include <vector>



#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"
#include "flint/arith.h"
#include "lib/near_entropic.h"
#include <stdio.h>
using std::cout, std::endl;



const int MAX_SYM = 16;
const int SEQ_LEN = 16;

// Compute the information content
#include <cmath>    // For std::log2 and std::round
double measureEntropy(const std::vector<uint8_t>& seq, int maxSym) {
	
    double n = static_cast<double>(seq.size());
    
    // Handle empty string case to avoid division by zero
    if (n == 0.0) return 0.0;

    // Pre-calculate frequencies (Optimization: O(N) instead of O(N^2))
    std::vector<uint8_t> counts(maxSym);	
    for (uint8_t val : seq) {    
        counts[val]++;
    }

    double total_log_sum = 0.0;
    for (uint8_t val : seq) {           
        // Probability of the current symbol
        double p = counts[val] / n;
        total_log_sum += std::log2(p);
    }

    // Apply negative sign and round to 2 decimal places
    // C++ standard round() rounds to nearest integer, so we multiply/divide by 100
    double result = -total_log_sum;
    return std::round(result * 100.0) / 100.0;
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
	fmpz_mat_print_pretty(mat);
}

int main() {
		
	//std::vector<uint8_t> valSeq = {1,7,7,1,14,7,0,11,2,13,13,12,11,2,0,7};
	//double ent = measureEntropy(valSeq, MAX_SYM);
	//cout << "ent: " << ent << endl;
	fmpz_mat_t table;
	gen_rgf_table(4, 2, table);
	serialize_mat("test.dat", table);
	fmpz_mat_clear(table);
	
	fmpz_mat_t table2;
	deserialize_mat("test.dat", table2);
	fmpz_mat_clear(table2);
	
	//fmpz_t entRank;
	//fmpz_init(entRank);
	//near_entropic_rank(valSeq, MAX_SYM, entRank);			
	//cout << "---" << endl;
	//std::vector<uint8_t> rgfSeq(SEQ_LEN);
	//near_entropic_unrank(entRank, SEQ_LEN, MAX_SYM, rgfSeq);	
	//printVector(rgfSeq);
	//fmpz_clear(entRank);
	return 0;
}