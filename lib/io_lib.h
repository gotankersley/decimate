#pragma once
#include <iostream>
#include <vector>
#include <cstdint>
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"


void printVector(std::vector<uint8_t>& vals);
void printVector(std::vector<int>& vals);

void serialize_mat(const char* filename, fmpz_mat_t mat);
void deserialize_mat(const char* filename, fmpz_mat_t mat); 