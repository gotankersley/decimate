#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include <stdio.h>
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"

void gen_rgf_table(int n, int k, fmpz_mat_t tableOut);
void serialize_mat(const char* filename, fmpz_mat_t mat);
void deserialize_mat(const char* filename, fmpz_mat_t mat); 

void rgf_rank(std::vector<uint8_t>& rgf, int k, fmpz_t rankOut);
void rgf_rank_table(std::vector<uint8_t>& rgf, int k, fmpz_mat_t table, fmpz_t rankOut);

void rgf_unrank(fmpz_t rank, int n, int k, std::vector<uint8_t>& rgfOut);
void rgf_unrank_table(fmpz_t rank, int n, int k, fmpz_mat_t table, std::vector<uint8_t>& rgfOut);