#include "nearer_entropic.h"
using std::cout, std::endl;


const bool DEBUG = true;
const int INVALID = -1;

//Combinatorial function to count the ways to rank by symbol count - defined in near_entropic.cpp
void addSymbolSection(fmpz_t rank, int seqLen, int maxSym, int k);

/*
#ABOUT:  This is a function to rank a sequence in sorted entropic order, rather than the regular place-value base system.
# This means that sequences with the less entropy, (i.e. fewer distinct symbols), like 'AAA', 'BBB', are ranked first,
# and those with the more entropy, (i.e. more differing symbols), like 'ABCD', are ranked last.
#
# Sequences are uniquely identified in the following way:
#  1. By symbol section, (i.e. all the sequences with only 1 symbol are grouped together, and so on)
#  2. By Set Partition largest part section			--|
#  3. By Set Partition initial part section			  +-- 2,3,4 Are in lieu of an RGF
#  4. By Set Partition initial element combination	--|
#  -------------	
#	Note: Everything after this has the exact same amount of entropy, 
#         so the rest is just one possible scheme to identify the sequence
#  -------------
#  5. By the combination rank of the symbols used from the set of total available symbols
#  6. By the permutation rank of the mapping the symbols to the set partition
*/
void nearer_entropic_rank(std::vector<uint8_t>& valSeq, int maxSym, fmpz_t rankOut) {

	//Pre-process seq	
	int seqLen = valSeq.size();
		
	std::vector<int> valToSym(maxSym); //Put val in get sym out	
	for (int i = 0; i < maxSym; i++) {
		valToSym[i] = INVALID;
	}

	std::vector<std::set<int>> setPart;	
	indexed_set_t unusedElements;
	
	int symCount = INVALID;	 
	for (int i = 0; i < seqLen; i++){
		uint8_t val = valSeq[i];
		uint8_t sym;
		//Create map of seq -> vals
		if (valToSym[val] == INVALID) { //First time this symbol has been seen		
			symCount++;
			valToSym[val] = symCount;
			sym = symCount;	
			std::set<int> part = {i};
			setPart.push_back(part);
		}
		else {
			sym = valToSym[val];
			setPart[sym].insert(setPart[sym].end(), i);
		}
		unusedElements.insert(i);
	}
	symCount++; //So that it reflects the total properly 
	
	if (DEBUG) {
		cout << "Ranking: ";		
		printVector(valSeq);
		printSetPart(setPart);
	}
	
	//Permutation of symbols -> vals
	std::vector<uint8_t> vals;
	std::vector<uint8_t> symPerm;	
	
	for (int i = 0; i < maxSym; i++) { 
		//Traverse vals in order
		if (valToSym[i] != INVALID) {
			vals.push_back(i);			
			symPerm.push_back(valToSym[i]);
		}		
	}
	
	if (DEBUG) {
		cout << "Sym Count: " << symCount << endl;		
		cout << "Vals: ";
		printVector(vals);		
	}
	
	
	fmpz_zero(rankOut);	
	
	// 1. Add symbol sections
	for (int i = 1; i < symCount; i++) {				
		addSymbolSection(rankOut, seqLen, maxSym, i);		
	}
		
	if (DEBUG) {
		cout << "Rank after symbol section: ";
		fmpz_print(rankOut);
		cout << endl << endl;
	}
		
		
	//NOTE: Ranking the set-partition here is the real heart of this algorithm.  
	//As the CAGES book says on ranking Set Partitions: "It is both convenient and natural to use an RGF for ranking".
	//Unfortunately, this method doesn't use an RGF, (and, thus, is neither convenient or natural).  
	//However, the extra complexity is required to achieve a better entropic order, as RGF order is not entropic
	int n = seqLen; 
	fmpz_t stirRank;
	fmpz_t count;
	fmpz_t elementSectionSize;
	fmpz_t elementRank;
	
	fmpz_init(stirRank);		
	fmpz_init(count);		
	fmpz_init(elementSectionSize);		
	fmpz_init(elementRank);		
	
	fmpz_zero(stirRank);
		
	int prevLargestPartSize = (n-symCount) + 1;	
	
	
	//Loop minus-1, because there is only a single way to do the last part that uses all the remaining elements
	for (int p = 0; p < symCount-1; p++) { 		
		int k = symCount - p;
		int m = get_size_of_largest_part(setPart, p);
		int r = setPart[p].size();
		if (DEBUG) {
			cout << "---- SET PART ITERATION: " << p << " ---" << endl;
			cout << "K: " << k << endl;
			cout << "M: " << m << endl;
			cout << "R: " << r << endl;
		}
		//std::string filename = "K" + std::to_string(k) + ".mat";		
		fmpz_mat_t coeffs;	
		gen_coeff_col(n, k, prevLargestPartSize, coeffs);
		
		//deserialize_mat(filename.c_str(), coeffs); //Load precalculated coefficients from file
		cout << "Coeffs: " << n << "," << k << "," << prevLargestPartSize << endl;
		fmpz_mat_print_pretty(coeffs);
		cout << endl << endl;
		
		// 2. By Set Partition largest part section	
		if (DEBUG) {
			cout << "Stir Rank before Largest part section: ";
			fmpz_print(stirRank);
			cout << endl;
		}
		for (int i = prevLargestPartSize; i > m; i--) {
			stirling2_max(coeffs, n, k, i, count);
			fmpz_add(stirRank, stirRank, count);
		}
		prevLargestPartSize = m;	
		if (DEBUG) {
			cout << "Stir Rank after Largest part section: ";
			fmpz_print(stirRank);
			cout << endl;
		}		
		//Note: At this point in the unranking we have just found the value of m (largest part size)
		
		// 3. By Set Partition initial part size section			
		for (int i = m; i > r; i--) {
			stirling2_max_r(coeffs, n, k, m, i, count);
			//Count may be zero, but it will just be ignored
			fmpz_add(stirRank, stirRank, count);			
		}
		//Note: At this point in the unranking we know the length of the initial part of the set partition r
		stirling2_max_r(coeffs, n, k, m, r, count); //Count here is initialPartSectionSize
		fmpz_mat_clear(coeffs);
		
		// 4. By Set Partition initial element combination 		
		unusedElements.erase(unusedElements.begin());		
		if (r > 1) { //If r == 1, then it only has a single element, and there's only one way
			//WARNING: This is a particularly tricky bit where we need to use an indexed position
			//of our set of elements, and then remove the elements, which changes all the upstream
			//index positions. (But, only AFTER they've all been used)
			//NOTE:  To solve this, we are using the __gnu_pbds extention which implements
			//an indexed ordered set: https://www.geeksforgeeks.org/cpp/ordered-set-gnu-c-pbds/			
			//TODO: Perhaps this could be optimized by using Factoradics instead?
			fmpz_bin_uiui(elementSectionSize, n-1, r-1); //The first element is always fixed as the lowest		
			fmpz_tdiv_q(elementSectionSize, count, elementSectionSize);
			
			std::vector<int> elementIds(r-1);
			int e = 0;			
			std::set<int>::iterator it = setPart[k-1].begin();						
			for (it++; it != setPart[k-1].end(); it++) { //Ignore first element					
				int elementVal = *it;
				int elementId = unusedElements.order_of_key(elementVal);
				elementIds[e] = elementId;
				//cout << "elementVal: " << +elementVal << endl;
				//cout << "elementId: " << +elementId << endl;
				//cout << "elementIds: " << elementIds[e] << endl;
				e++;
			}			
			
			//Remove all the used elements from the set
			it = setPart[k-1].begin();
			for (it++; it != setPart[k-1].end(); it++) { //Ignore first element			
				int elementVal = *it;
				unusedElements.erase(unusedElements.find(elementVal));							
			}
			
			if (DEBUG) {
				cout << "Element Ids: ";
				printVector(elementIds);				
				cout << endl << endl;
			}
			
			comb_rank(elementIds, elementRank);
			fmpz_add(stirRank, stirRank, elementRank);			
		}
		//Iterate - essentially we are recursing on all the rest of the parts
		n -= r;
	}
	fmpz_clear(count);
	fmpz_clear(elementSectionSize);
	fmpz_clear(elementRank);
	
	//Calculate section sizes	
	fmpz_t stirSectionSize;
	fmpz_init(stirSectionSize);
	uint64_t totalComb = comb(maxSym, symCount);
	fmpz_t combSectionSize;
	fmpz_init(combSectionSize);
	fmpz_fac_ui(combSectionSize, symCount);
	fmpz_mul_ui(stirSectionSize, combSectionSize, totalComb);
	
	if (DEBUG) {
		cout << "Comb Section Size: ";
		fmpz_print(combSectionSize);
		cout << endl;
		cout << "Stir Section Size: ";
		fmpz_print(stirSectionSize);
		cout << endl;
	}
	fmpz_addmul(rankOut, stirRank, stirSectionSize);
	fmpz_clear(stirSectionSize);
	
	if (DEBUG) {
		cout << "Stir Rank: ";
		fmpz_print(stirRank);
		cout << endl << endl;
	}
	
	fmpz_clear(stirRank);
	
	
	//  5. Add the combination rank 
	uint64_t combRank = comb_rank(vals);	
	fmpz_addmul_ui(rankOut, combSectionSize, combRank);
	if (DEBUG) {
		cout << "Comb Rank: " << combRank << endl;		
		cout << "Comb Vals: ";
		printVector(vals);
		cout << "Comb Section Size: ";
		fmpz_print(combSectionSize);
		cout << endl << endl;
	}
	fmpz_clear(combSectionSize);
	
	//  6. Add the Sym Perm Rank (Myrvold)	
	fmpz_t symRank;
	fmpz_init(symRank);
	if (DEBUG) {
		cout << "Sym Perm: ";
		printVector(symPerm);
	}
	myrvold_rank(symPerm, symRank);
	fmpz_add(rankOut, rankOut, symRank);
	
	if (DEBUG) {			
		cout << "Sym Rank: ";
		fmpz_print (symRank);
		cout << endl;
		
		cout << "Final Rank: " << endl;
		fmpz_print(rankOut);
		cout << endl;
	}
	fmpz_clear(symRank);	
	
}

void nearer_entropic_unrank(fmpz_t rank, int seqLen, int maxSym, std::vector<int>& countsOut, std::vector<uint8_t>& valSeqOut) {
			
	
	//  1. Get symbol section
	if (DEBUG) {
		cout << "Unranking" << endl;
		cout << "Rank: " << endl;
		fmpz_print(rank);
		cout << endl;
	}
	
	fmpz_t count;
	fmpz_t prevCount;
	
	fmpz_init(count);
	fmpz_init(prevCount);
	
	fmpz_zero(count);
	fmpz_zero(prevCount);
	int symCount = INVALID;
	for (int i = 1; i <= maxSym; i++) {
		addSymbolSection(count, seqLen, maxSym, i);		
		if (fmpz_cmp(count, rank) > 0) {		
			symCount = i;
			fmpz_sub(rank, rank, prevCount);
			break;
		}
		else fmpz_set(prevCount, count);
	}
	
	if (DEBUG) {
		cout << "Sym Count: " << symCount << endl;		
	}	
	
	
	//Calculate section sizes	
	uint64_t totalComb = comb(maxSym, symCount);
	fmpz_t combSectionSize;
	fmpz_init(combSectionSize);
	fmpz_fac_ui(combSectionSize, symCount);
	
	fmpz_t stirSectionSize;
	fmpz_init(stirSectionSize);
	fmpz_mul_ui(stirSectionSize, combSectionSize, totalComb);
	

	if (DEBUG) {
		cout << "Comb Section Size: " << endl;
		fmpz_print(combSectionSize);
		cout << endl;
		cout << "Stir Section Size: " << endl;
		fmpz_print(stirSectionSize);
		cout << endl;
	}
	
	// 2. Get the values from the combination rank of symbols
	fmpz_t rankModStir;	
	fmpz_init(rankModStir);	
	fmpz_mod(rankModStir, rank, stirSectionSize);			
	fmpz_tdiv_q(rankModStir, rankModStir, combSectionSize);	
	uint64_t combRank = fmpz_get_ui(rankModStir);	
	fmpz_clear(rankModStir);	
	std::vector<uint8_t> combVals(symCount);	
	comb_unrank(combRank, maxSym, symCount, combVals);
	
	
	if (DEBUG) {
		cout << "Comb Rank: " << combRank << endl;
		cout << "Comb Vals: ";
		printVector(combVals);
	}
	
	//  3. Get the Sym Perm from Myrvold Rank	
	fmpz_t symRank;
	fmpz_init(symRank);
	fmpz_mod(symRank, rank, combSectionSize);	
	std::vector<uint8_t> symPerm(symCount);
	if (DEBUG) {	
		cout << "Sym Rank: ";
		fmpz_print (symRank);
		cout << endl;
	}
	myrvold_unrank(symRank, symPerm);
	fmpz_clear(symRank);
	fmpz_clear(combSectionSize);
	std::vector<uint8_t> invPerm(symCount);
	for (int i = 0; i < symCount; i++) {
		invPerm[symPerm[i]] = i;
	}
	
	if (DEBUG) {
		cout << "Sym Perm: ";
		printVector(symPerm);		
	}

	
	// 4. Get the Set Partition from Stirling2 rank		
	// Note: We're also applying inverse perm to recreate seq	
	fmpz_t stirRank;
	fmpz_init(stirRank);
	fmpz_tdiv_q(stirRank, rank, stirSectionSize);
	fmpz_clear(stirSectionSize);	

	if (DEBUG) {
		cout << "Stir Rank: " << endl;
		fmpz_print(stirRank);
		cout << endl;					
	}	
	// Unrank the Set-Partition one part at a time
	int n = seqLen;
	int prevLargestPartSize = (n-symCount)+1;
	indexed_set_t unusedElements;
	for (int i = 0; i < n; i++) {
		unusedElements.insert(i);
	}
	fmpz_t count2;
	fmpz_t elementSectionSize;
	fmpz_t elementRank;
	
	fmpz_init(count2);	
	fmpz_init(elementSectionSize);	
	fmpz_init(elementRank);	
	
	for (int p = 0; p < symCount-1; p++) {
		int k = symCount - p;
		std::string filename = std::to_string(k) + ".mat";		
		fmpz_mat_t coeffs;
		deserialize_mat(filename.c_str(), coeffs); //Load precalculated coefficients from file
				
		uint8_t partVal = combVals[invPerm[k]];		
				
		//Set Partition largest part section	
		fmpz_zero(count);
		fmpz_zero(prevCount);
		int m;
		for (m = prevLargestPartSize; m > 0; m--) {
			stirling2_max(coeffs, n, k, m, count2);
			fmpz_add(count, count, count2);
			if (fmpz_cmp(count, stirRank) > 0) {
				fmpz_sub(stirRank, stirRank, prevCount);
				break;
			}
			else fmpz_set(prevCount, count);			
		}		
		prevLargestPartSize = m;
		
		//Set Partition initial part size section	
		fmpz_zero(count);
		fmpz_zero(prevCount);
		int r;
		int lastR;
		for (r = m; r > 0; r--) {
			stirling2_max_r(coeffs, n, k, m, r, count2);
			lastR = r;
			//Count may be zero, but it will just be ignored
			fmpz_add(count, count, count2);
			
			if (fmpz_cmp(count, stirRank) > 0) {
				fmpz_sub(stirRank, stirRank, prevCount);
				break;
			}
			else fmpz_set(prevCount, count);
		}
		if (lastR != r) stirling2_max_r(coeffs, n, k, m, r, count2); //initialPartSectionSize
		fmpz_mat_clear(coeffs);
		countsOut[k] = r;
		
		//Set Partition initial element combination.
		//The first element is always fixed as the lowest
		fmpz_bin_uiui(elementSectionSize, n-1, r-1); //The first element is always fixed as the lowest		
		fmpz_tdiv_q(elementSectionSize, count2, elementSectionSize);
		fmpz_tdiv_q(elementRank, stirRank, elementSectionSize);
		
		std::vector<int> elementIds(r-1);	
		comb_unrank(elementRank, n-1, r-1, elementIds);
		int initialElement = *unusedElements.begin();
		unusedElements.erase(unusedElements.begin());		
		valSeqOut[initialElement] = partVal;
		if (r > 1) {
			for (size_t i = 0; i < elementIds.size(); i++) {
				int elementId = elementIds[i];
				int element = *unusedElements.find_by_order(elementId);
				
				//Apply values here so we don't have to do another loop over the sequence				
				valSeqOut[element] = partVal;								
			}
			//Todo store iterators from previous? std::vector<indexed_set_t::iterator> its
			for (size_t i = 0; i < elementIds.size(); i++) {
				int elementId = elementIds[i];
				unusedElements.erase(unusedElements.find_by_order(elementId));				
			}
		}
		fmpz_mul(elementRank, elementSectionSize, elementRank);
		fmpz_sub(stirRank, stirRank, elementRank);
		n -= r; //Iterate		
	}	
	fmpz_clear(count);
	fmpz_clear(count2);
	fmpz_clear(prevCount);
	fmpz_clear(elementSectionSize);
	fmpz_clear(elementRank);
	fmpz_clear(stirRank);
	
}



