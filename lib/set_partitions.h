#pragma once
#include <vector>
#include <cstdint>
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"

void gen_rgf_table(int n, int k, fmpz_mat_t tableOut);

void rgf_rank(std::vector<uint8_t>& rgf, int k, fmpz_t rankOut);

void rgf_unrank(fmpz_t rank, int n, int k, std::vector<uint8_t>& rgfOut);