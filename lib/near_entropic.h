#pragma once
#include <iostream>
#include <cstdint>
#include <cmath>
#include <vector>

#include "flint/fmpz.h"
#include "flint/arith.h"
#include "combinations.h"
#include "permutations.h"
#include "rgf.h"


void printVector(std::vector<uint8_t>& vals);
double measureEntropy(std::vector<int>& counts, int seqLen);

void near_entropic_rank(std::vector<uint8_t>& valSeq, int maxSym, fmpz_t rankOut);
void near_entropic_unrank(fmpz_t rank, int seqLen, int maxSym, std::vector<int>& countsOut, std::vector<uint8_t>& rgfOut);