#pragma once
#include <vector>
#include <cstdint>
#include "flint/fmpz.h"
#include "flint/arith.h"

uint64_t comb(int n, int k); //Binomial co-efficient of n-choose-k	

uint64_t comb_rank(std::vector<uint8_t>& vals);
void comb_rank(std::vector<int>& vals, fmpz_t rankOut);

void comb_unrank(uint64_t rank, int n, int k, std::vector<uint8_t>& valsOut);
void comb_unrank(fmpz_t rank, int n, int k, std::vector<int>& valsOut);
	