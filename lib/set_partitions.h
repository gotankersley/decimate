#pragma once
#include <vector>
#include <string>
#include <set>
#include <cstdint>
#include <stdio.h>
#include "io_lib.h"
#include "combinations.h"
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"
#include "flint/arith.h"
#include "flint/fmpq.h"
#include "flint/fmpq_poly.h"


void gen_initial_coeffs(int m, fmpz_mat_t coeffsOut);
void gen_coeff_table(int maxN, int maxK, int m);//Output is serialized to files

//TODO - cache coeffs
void stirling2_max_lt(int n, int k, int m, fmpz_t countOut); //Polynomial EGF
void stirling2_max_between(int n, int k, int mHi, int mLo, fmpz_t countOut);


void stirling2_max_initial_lt(int n, int k, int m, int r, fmpz_t countOut); //Q-nomial EGF
void stirling2_max_initial_ge(int n, int k, int m, int r, fmpz_t countOut);
void stirling2_max_initial_gt(int n, int k, int m, int r, fmpz_t countOut);

void printSetPart(std::vector<std::set<int>>& setPart);
int get_size_of_largest_part(std::vector<std::set<int>>& setPart, int start);
