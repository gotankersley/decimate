#pragma once
#include <vector>
#include <string>
#include <set>
#include <cstdint>
#include <stdio.h>
#include "io_lib.h"
#include "combinations.h"
#include "flint/arith.h"
#include "flint/fmpz.h"
#include "flint/fmpq.h"
#include "flint/fmpz_mat.h"
#include "flint/fmpq_mat.h"
#include "flint/fmpq_poly.h"



void gen_k_facts(int k, fmpz_mat_t kFactsOut); 
void gen_coeffs(int m, fmpq_mat_t coeffsOut); 



void stirling2_max_lt (int n, int k, int m, fmpz_t nFact, fmpz_t kFact, fmpq_mat_t coeffs, fmpz_t countOut); //Polynomial EGF
void stirling2_max_between(int n, int k, int mHi, int mLo, fmpz_t nFact, fmpz_t kFact, fmpq_mat_t coeffs, fmpz_t countOut);


void stirling2_max_initial_lt(int n, int k, int m, int r, fmpz_t nFact, fmpz_t kFact, fmpq_mat_t coeffs, fmpz_t countOut); //Q-nomial EGF
void stirling2_max_initial_ge(int n, int k, int m, int r, fmpz_t nFact, fmpz_t kFact, fmpq_mat_t coeffs, fmpz_t countOut);
void stirling2_max_initial_gt(int n, int k, int m, int r, fmpz_t nFact, fmpz_t kFact, fmpq_mat_t coeffs, fmpz_t countOut);

void printSetPart(std::vector<std::set<int>>& setPart);
int get_size_of_largest_part(std::vector<std::set<int>>& setPart, int start);

/*
//Old / Non-Precached
void stirling2_max_initial_lt_old(int n, int k, int m, int r, fmpz_t nFact, fmpz_t kFact, fmpz_t countOut); //Q-nomial EGF
void stirling2_max_lt_old(int n, int k, int m, fmpz_t nFact, fmpz_t kFact, fmpz_t countOut); //Polynomial EGF
*/