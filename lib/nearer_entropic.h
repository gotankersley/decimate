#pragma once
#include <iostream>
#include <cstdint>
#include <vector>
#include <set>
#include <ext/pb_ds/assoc_container.hpp> //Indexed set
#include <ext/pb_ds/tree_policy.hpp> //Indexed set

#include "flint/fmpz.h"
#include "flint/arith.h"
#include "base_lib.h"
#include "io_lib.h"
#include "combinations.h"
#include "permutations.h"
#include "set_partitions.h"

using namespace __gnu_pbds; //Indexed Set
#define indexed_set_t tree<int, null_type,std::less<int>, rb_tree_tag,tree_order_statistics_node_update>

void nearer_entropic_rank(std::vector<uint8_t>& valSeq, int maxSym, fmpz_t rankOut);
void nearer_entropic_unrank(fmpz_t rank, int seqLen, int maxSym, std::vector<int>& countsOut, std::vector<uint8_t>& rgfOut);