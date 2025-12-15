#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include <stdio.h>
#include "combinations.h"
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"
#include "flint/arith.h"


const int NUM = 0;
const int DEN = 1;
const int FRACTION_SIZE = 2;

void mult_poly(fmpz_mat_t coeffs1, fmpz_mat_t coeffs2, fmpz_mat_t coeffsOut);
void gen_initial_coeffs(int m, fmpz_mat_t coeffsOut);
void gen_coeff_table(int maxN, int maxK, int m);//Output is serialized to files

void serialize_mat(const char* filename, fmpz_mat_t mat);
//void deserialize_mat(const char* filename, fmpz_mat_t mat); 


//seq_to_set_part

void stirling2_max_less_than_table(fmpz_mat_t coeffs, int n, int k, int m, fmpz_t countOut);
void stirling2_max(fmpz_mat_t coeffs, int n, int k, int m, fmpz_t countOut);
void stirling2_max_r(fmpz_mat_t coeffs, int n, int k, int m, int r, fmpz_t countOut);




