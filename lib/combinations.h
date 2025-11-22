#pragma once
#include <vector>
#include <cstdint>

uint64_t comb(int n, int k); //Binomial co-efficient of n-choose-k	

uint64_t comb_rank(std::vector<uint8_t> vals);
std::vector<uint8_t> comb_unrank(uint64_t rank, int n, int k);
	