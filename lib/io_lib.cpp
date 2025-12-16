#include "io_lib.h"

void printVector(std::vector<uint8_t>& vals) {    
for (int val : vals) {
		std::cout << +val << ","; 
	}
	std::cout << std::endl;
}

void printVector(std::vector<int>& vals) {    
for (int val : vals) {
		std::cout << +val << ","; 
	}
	std::cout << std::endl;
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