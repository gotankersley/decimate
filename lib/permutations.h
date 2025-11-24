#pragma once
#include <vector>
#include <cstdint>
#include "flint/fmpz.h"

void myrvold_rank_recur(int n, std::vector<uint8_t>& perm, std::vector<uint8_t>& invPerm, fmpz_t rankOut);
void myrvold_rank(std::vector<uint8_t> perm, fmpz_t rankOut);

void myrvold_unrank_recur(fmpz_t rank, int n, std::vector<uint8_t>& permOut);
void myrvold_unrank(fmpz_t rank, std::vector<uint8_t>& permOut);

