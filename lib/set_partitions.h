#pragma once
#include <vector>
#include <cstdint>
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"

void gen_rgf_table(fmpz_mat_t table, int n, int k);

void rgf_rank(fmpz_t rank, std::vector<uint8_t>& rgf, int k);

void rgf_unrank(std::vector<uint8_t>& rgf, fmpz_t rank, int n, int k);