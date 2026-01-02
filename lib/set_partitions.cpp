#include "set_partitions.h"
#include <iostream>

using std::cout, std::endl;

void gen_k_facts(int k, fmpz_mat_t kFactsOut) {
	fmpz_mat_init(kFactsOut, 1, k+1);
	
	for (int i = 1; i <= k; i++) {			
		fmpz_fac_ui(fmpz_mat_entry(kFactsOut, 0, i), i);		
	}
}



void stirling2_max_lt(int n, int k, int m, fmpz_t kFact, fmpz_t countOut) { //Polynomial EGF
	fmpz_zero(countOut);
	if (k > n) return;
	else if (m <= 0) return;
	//TODO - cache coeffs
	
	// 1. Create coefficients for polynomial B(x) = x/1! + x^2/2! + ... + x^m/m!
	fmpz_t coeffNum;
	fmpz_t coeffDen;	
	
	fmpq_t coeff;
	fmpq_poly_t poly;
	
	fmpz_init(coeffNum);
	fmpz_init(coeffDen);	
	fmpq_init(coeff);
	fmpq_poly_init(poly); //TODO: fmpq_poly_init2 instead?
	
	fmpz_one(coeffNum);
	fmpz_one(coeffDen);
	for (int i = 1; i < m+1; i++) {
		fmpz_mul_ui(coeffDen, coeffDen, i);
		fmpq_set_fmpz_frac(coeff, coeffNum, coeffDen);		
		fmpq_poly_set_coeff_fmpq(poly, i-1, coeff);
	}
	
	
	// 2. Polynomial power to perform convolution
    // NOTE: We factored out x^k, so coefficient for x^n is at index n - k
	int targetIndex = n - k;
	
	fmpq_poly_pow_trunc(poly, poly, k, targetIndex+1);
	
	
	// 3. Extract coefficient
    if (targetIndex >= 0 && targetIndex < fmpq_poly_length(poly)) {
		fmpq_poly_get_coeff_fmpq(coeff, poly, targetIndex);		
		
		// 4. Multiply by n! / k!
		//TODO: Precalc
		fmpz_t factorialN;			
		
		fmpz_init(factorialN);		
		
		fmpz_fac_ui(factorialN, n);			
		
		fmpz* numPtr = fmpq_numref(coeff);
		fmpz* denPtr = fmpq_denref(coeff);
		
		fmpz_mul(coeffNum, factorialN, numPtr);		
		fmpz_mul(coeffDen, kFact, denPtr);
		fmpz_tdiv_q(countOut, coeffNum, coeffDen);
		
		fmpz_clear(factorialN);		
	}
	
	//Clear
	fmpq_clear(coeff);
	fmpz_clear(coeffNum);
	fmpz_clear(coeffDen);
	fmpq_poly_clear(poly);
}
void stirling2_max_between(int n, int k, int mHi, int mLo, fmpz_t kFact, fmpz_t countOut) {
	//This calculates the Stirling2 number of set-parts with n elements, k parts, 
	//and having a max part size that is between mHi and mLo 
	fmpz_t loCount;
	fmpz_init(loCount);
	stirling2_max_lt(n, k, mHi, kFact, countOut);
	stirling2_max_lt(n, k, mLo, kFact, loCount);	
	fmpz_sub(countOut, countOut, loCount);
	fmpz_clear(loCount);
}


void stirling2_max_initial_lt(int n, int k, int m, int r, fmpz* kFact, fmpz_t countOut) { //Q-nomial EGF
	//This calculates the Stirling2 number of set-parts with n elements, k parts, 
	//and having a max part size that is exactly m, and that has an initial part size less-than r
	fmpz_zero(countOut);
	// Edge Cases
	if (k < 0) return;
	else if (k > n) return;
	else if (n == 0) {
		if (k == 0) fmpz_one(countOut);
		return;		
	}
	else if (k == 0) return;
	else if (r > m) return;
	
	// Create Polynomials
	fmpz_t coeffNum;
	fmpz_t coeffDen;	
	fmpz_t prevDen;	
	fmpq_t coeff;
	fmpq_poly_t poly;
	fmpq_poly_t polyQ;
	
	fmpz_init(coeffNum);
	fmpz_init(coeffDen);	
	fmpz_init(prevDen);	
	fmpq_init(coeff);
	fmpq_poly_init(poly); 
	fmpq_poly_init(polyQ); 
	
	fmpz_one(coeffNum);
	fmpz_one(coeffDen);
	fmpz_one(prevDen);
	for (int i = 0; i < m; i++) {
		fmpz_mul_ui(coeffDen, coeffDen, i+1);
		fmpq_set_fmpz_frac(coeff, coeffNum, coeffDen);
		fmpq_poly_set_coeff_fmpq(poly, i, coeff);
		
		// Check if this power is valid
		if (i >= (r-1)) {
			fmpq_set_fmpz_frac(coeff, coeffNum, prevDen);
			fmpq_poly_set_coeff_fmpq(polyQ, i, coeff);
		}
		fmpz_set(prevDen, coeffDen);
	}
	
	// Compute powers
	int targetIndex = n - k;
	if (k == 1) {
		fmpq_poly_one(poly);	
	}
	else {
		fmpq_poly_pow_trunc(poly, poly, k-1, targetIndex+1); //TODO: Should this be truncated?
	}
	//Multiply polynomials
	fmpq_poly_mullow(polyQ, polyQ, poly, targetIndex+1);
	
	//Extract coefficient
    if (targetIndex >= 0 && targetIndex < fmpq_poly_length(polyQ)) {
		fmpq_poly_get_coeff_fmpq(coeff, polyQ, targetIndex);
		
		// 4. Multiply by (n-1)! / (k-1)!
		//TODO: Precalc factorials
		fmpz_t factorialN;	
		//fmpz_t factorialK;	
		
		fmpz_init(factorialN);
		//fmpz_init(factorialK);
		
		fmpz_fac_ui(factorialN, n-1);	
		//fmpz_fac_ui(factorialK, k-1);	
		
		fmpz* numPtr = fmpq_numref(coeff);
		fmpz* denPtr = fmpq_denref(coeff);
		
		fmpz_mul(coeffNum, factorialN, numPtr);
		//fmpz_mul(coeffDen, factorialK, denPtr);
		fmpz_mul(coeffDen, kFact, denPtr);
		fmpz_tdiv_q(countOut, coeffNum, coeffDen);
		fmpz_clear(factorialN);
		//fmpz_clear(factorialK);		
	}
	
	
	//Clear
	fmpq_clear(coeff);
	fmpz_clear(coeffNum);
	fmpz_clear(coeffDen);
	fmpz_clear(prevDen);
	fmpq_poly_clear(poly);
	fmpq_poly_clear(polyQ);
	
}
void stirling2_max_initial_ge(int n, int k, int m, int r, fmpz* kFact, fmpz_t countOut) {
	//This calculates the Stirling2 number of set-parts with n elements, k parts, 
	//and having a max part size that is exactly m, and that has an initial part size greater-than-or-equal to r
	fmpz_t mCount;
	fmpz_init(mCount);
	stirling2_max_initial_lt(n, k, m, r, kFact, countOut);
	stirling2_max_initial_lt(n, k, m-1, r, kFact, mCount);	
	fmpz_sub(countOut, countOut, mCount);
	fmpz_clear(mCount);
}
void stirling2_max_initial_gt(int n, int k, int m, int r, fmpz* kFact, fmpz_t countOut) {
	//This calculates the Stirling2 number of set-parts with n elements, k parts, 
	//and having a max part size that is exactly m, and that has an initial part size greater-than r
	fmpz_t mCount;
	fmpz_init(mCount);
	stirling2_max_initial_lt(n, k, m, r+1, kFact, countOut);
	stirling2_max_initial_lt(n, k, m-1, r+1, kFact, mCount);	
	fmpz_sub(countOut, countOut, mCount);
	fmpz_clear(mCount);
}

/*

void gen_initial_coeffs(int m, fmpz_mat_t coeffsOut) {
	//Create coefficients for polynomial B(x) = x/1! + x^2/2! + ... + x^m/m!
	//Table dimensions: 
	// - two rows, first is the numerator, second is the denominator
	// - column index i corresponds to power x^i	
	fmpz_mat_init(coeffsOut, FRACTION_SIZE, m);	
	
	fmpz_t fact;
	fmpz_init(fact);
	fmpz_one(fact);
	for (int i = 0; i < m; i++) {
		fmpz_one(fmpz_mat_entry(coeffsOut, NUM, i));
		fmpz_mul_ui(fact, fact, i+1);
		fmpz_set(fmpz_mat_entry(coeffsOut, DEN, i), fact);		
	}
	fmpz_clear(fact);
}


void gen_coeff_table(int maxN, int maxK, int m) {//Output is serialized
	//Note:  This is the precalculated part that is used by the
	//stirling2_max_less_than() function to make it viable.
	//We are calculating everything but the final step where n! is added
	fmpz_mat_t baseCoeffs;	
	gen_initial_coeffs(m, baseCoeffs);	
	serialize_mat("K1.mat", baseCoeffs);	
	
	fmpz_mat_t curCoeffs;		
	fmpz_mat_init(curCoeffs, FRACTION_SIZE, m);
	fmpz_mat_set(curCoeffs, baseCoeffs); //Copy - used for successive powers
		
	fmpz_mat_t resCoeffs;	
	
	//Raise the base polynomial to the kth power
	fmpz_t factorialK;	
	fmpz_t gcd;
	
	fmpz_init(factorialK);
	fmpz_init(gcd);
	for (int k = 2; k <= maxK; k++) {
		fmpz_fac_ui(factorialK, k);			
		mult_poly(baseCoeffs, curCoeffs, resCoeffs);
		
		//Extract coefficients, and store - Note: we only need up to power of N-K
		int maxIdx = std::min(maxN-k+1, (int)(resCoeffs->c)); 
		fmpz_mat_t col;
		fmpz_mat_init(col, FRACTION_SIZE, maxIdx);
		for (int i = 0; i < maxIdx; i++) {
			//Simplify
			fmpz* n = fmpz_mat_entry(resCoeffs, NUM, i);
			fmpz* d = fmpz_mat_entry(resCoeffs, DEN, i);
			fmpz* colN = fmpz_mat_entry(col, NUM, i);
			fmpz* colD = fmpz_mat_entry(col, DEN, i);
						
			fmpz_mul(colD, d, factorialK); //Precalculate k!, so we don't have to do this later	
			
			fmpz_gcd(gcd, n, colD);
			if (fmpz_is_one(gcd)) {
				fmpz_set(colN, n);
			}
			else {
				fmpz_tdiv_q(colN, n, gcd);
				fmpz_tdiv_q(colD, colD, gcd);
			}			
			
		}
		
		std::string filename = "K" + std::to_string(k) + ".mat";
		serialize_mat(filename.c_str(), col);		
		fmpz_mat_clear(col);
		
		//Swap current
		fmpz_mat_swap(curCoeffs, resCoeffs);
		fmpz_mat_clear(resCoeffs);				
	}
	fmpz_clear(factorialK);
	fmpz_clear(gcd);
	
	
	fmpz_mat_clear(baseCoeffs);
	fmpz_mat_clear(curCoeffs);
	
}
void stirling2_max_less_than_coeffs(fmpz_mat_t coeffs, int n, int k, int m, fmpz_t countOut) {
	fmpz_init(countOut);
	fmpz_zero(countOut);
	if (k > n || m <= 0) return; //Count of 0
	
	int targetIdx = n-k;
	if (targetIdx < 0 || targetIdx >= coeffs->c) return; //Count of 0
		
	fmpz_t factorialN;
	fmpz_init(factorialN);
	fmpz_fac_ui(factorialN, n);
	
	fmpz* vn = fmpz_mat_entry(coeffs, NUM, targetIdx);
	fmpz* vd = fmpz_mat_entry(coeffs, DEN, targetIdx);
	
	fmpz_mul(countOut, factorialN, vn);
	fmpz_tdiv_q(countOut, countOut, vd); //Since this is a count, it will always be an integer
	fmpz_clear(factorialN);
	
}


void stirling2_max_eq(fmpz_mat_t coeffs, int n, int k, int m, fmpz_t countOut) {
	//This calculates the Stirling2 number of set-parts with n elements, k parts, and having a max part size that is exactly m
	fmpz_t count2;
	fmpz_init(count2);
	stirling2_max_less_than_coeffs(coeffs, n, k, m-1, count2);	
	stirling2_max_less_than_coeffs(coeffs, n, k, m, countOut);
	fmpz_sub(countOut, countOut, count2);
	fmpz_clear(count2);
}
*/



void printSetPart(std::vector<std::set<int>>& setPart) {
	cout << "Set Part: " << endl;
	for (int p = 0; p < (int)setPart.size(); p++) {
		cout << "Part: " << p << " [";
		std::set<int> part = setPart[p];
		
		for (const int element : part) {
			cout << +element << ",";         
		}		
		cout << "]" << endl;
	}
	std::cout << std::endl;
	cout << endl;
}

int get_size_of_largest_part(std::vector<std::set<int>>& setPart, int start) {
	int largestPart = -1;	
	
	for (int i = start; i < (int)setPart.size(); i++) {
		std::set<int> part = setPart[i];
		if ((int)part.size() > largestPart) largestPart = part.size();
	}	
	return largestPart;
}