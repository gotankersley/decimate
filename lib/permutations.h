#pragma once
#include <vector>
#include <cstdint>
#include "flint/fmpz.h"

void myrvold_rank_recur(fmpz_t rank, int n, std::vector<uint8_t>& perm, std::vector<uint8_t>& invPerm);
void myrvold_rank(fmpz_t rank, std::vector<uint8_t>& perm);

void myrvold_unrank_recur(fmpz_t rank, int n, std::vector<uint8_t>& perm);
void myrvold_unrank(std::vector<uint8_t>&perm, fmpz_t rank);

