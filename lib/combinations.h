#pragma once
#include <vector>
#include <cstdint>

uint64_t comb(int n, int k); //Binomial co-efficient of n-choose-k	

uint64_t comb_rank(std::vector<uint8_t>& vals);
void comb_unrank(std::vector<uint8_t>& vals, uint64_t rank, int n, int k);
	