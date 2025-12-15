#pragma once
#include <vector>
#include <cstdint>
#include <cmath>
#include "flint/fmpz.h"

void b2n(std::vector<uint8_t>& digits, int base, std::vector<int>& countsOut, fmpz_t nOut);

double measureEntropy(std::vector<int>& counts, int seqLen);