#include "set_partitions.h"
#include <iostream>

using std::cout, std::endl;

void gen_k_facts(int k, fmpz_mat_t kFactsOut) {
	fmpz_mat_init(kFactsOut, 1, k+1);
	
	for (int i = 1; i <= k; i++) {			
		fmpz_fac_ui(fmpz_mat_entry(kFactsOut, 0, i), i);		
	}
}

void gen_coeffs(int m, fmpq_mat_t coeffsOut) {
	fmpq_mat_init(coeffsOut, 1, m+1);	
	fmpz_t num;
	fmpz_t den;
	
	fmpz_init(num);
	fmpz_init(den);
	
	fmpz_one(num);
	fmpz_one(den);
	for (int i = 0; i <= m; i++) {		
		fmpz_mul_ui(den, den, i+1);		
		fmpq_set_fmpz_frac(fmpq_mat_entry(coeffsOut, 0, i), num, den);				
	}
	fmpz_clear(num);
	fmpz_clear(den);

}

void stirling2_max_lt(int n, int k, int m, fmpz_t nFact, fmpz_t kFact, fmpq_mat_t coeffs, fmpz_t countOut) { //Polynomial EGF
	fmpz_zero(countOut);
	if (k > n) return;
	else if (m <= 0) return;
	
	
	// Create coefficients for polynomial B(x) = x/1! + x^2/2! + ... + x^m/m!
	fmpq_poly_t poly;
	fmpq_poly_init(poly); //TODO: fmpq_poly_init2 instead?
	

	for (int i = 0; i < m; i++) {
		fmpq_poly_set_coeff_fmpq(poly, i, fmpq_mat_entry(coeffs, 0, i));
	}
	
	// Polynomial power to perform convolution
    // NOTE: We factored out x^k, so coefficient for x^n is at index n - k
	int targetIndex = n - k;
	
	fmpq_poly_pow_trunc(poly, poly, k, targetIndex+1);
	
	
	// Extract coefficient
    if (targetIndex >= 0 && targetIndex < fmpq_poly_length(poly)) {
		fmpq_t coeff;
		fmpz_t coeffNum;
		fmpz_t coeffDen;		
		
		fmpq_init(coeff);
		fmpz_init(coeffNum);
		fmpz_init(coeffDen);
		fmpq_poly_get_coeff_fmpq(coeff, poly, targetIndex);		
		
		// Multiply by n! / k!								
		fmpz* numPtr = fmpq_numref(coeff);
		fmpz* denPtr = fmpq_denref(coeff);
			
		
		fmpz_mul(coeffNum, nFact, numPtr);		
		fmpz_mul(coeffDen, kFact, denPtr);
		fmpz_tdiv_q(countOut, coeffNum, coeffDen);	
		
		fmpq_clear(coeff);
		fmpz_clear(coeffNum);
		fmpz_clear(coeffDen);		
	}
	
	//Clear
	fmpq_poly_clear(poly);
}
void stirling2_max_between(int n, int k, int mHi, int mLo, fmpz_t nFact, fmpz_t kFact, fmpq_mat_t coeffs, fmpz_t countOut) {
	//This calculates the Stirling2 number of set-parts with n elements, k parts, 
	//and having a max part size that is between mHi and mLo 
	fmpz_t loCount;
	fmpz_init(loCount);
	stirling2_max_lt(n, k, mHi, nFact, kFact, coeffs, countOut);
	stirling2_max_lt(n, k, mLo, nFact, kFact, coeffs, loCount);	
	fmpz_sub(countOut, countOut, loCount);
	fmpz_clear(loCount);
}

/*
void stirling2_max_initial_lt(int n, int k, int m, int r, fmpz_t nFact, fmpz_t kFact, fmpq_mat_t coeffs, fmpz_t countOut) { //Q-nomial EGF
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
	fmpq_poly_t poly;
	fmpq_poly_t polyQ;

	fmpq_poly_init(poly); 
	fmpq_poly_init(polyQ); 
	for (int i = 0; i < m; i++) {
		fmpq_poly_set_coeff_fmpq(poly, i, fmpq_mat_entry(coeffs, 0, i));
		// Check if this power is valid
		if (i >= (r-1)) {
			fmpq_poly_set_coeff_fmpq(polyQ, i, fmpq_mat_entry(coeffs, 0, i-1));
		}
	}
	
	// Compute powers
	int targetIndex = n - k;
	if (k == 1) {
		fmpq_poly_one(poly);	
	}
	else {
		fmpq_poly_pow_trunc(poly, poly, k-1, targetIndex+1); 
	}
	//Multiply polynomials
	fmpq_poly_mullow(polyQ, polyQ, poly, targetIndex+1);
	
	//Extract coefficient
    if (targetIndex >= 0 && targetIndex < fmpq_poly_length(polyQ)) {
		fmpz_t coeffNum;
		fmpz_t coeffDen;		
		fmpq_t coeff;
		
		fmpz_init(coeffNum);
		fmpz_init(coeffDen);		
		fmpq_init(coeff);
		
		fmpq_poly_get_coeff_fmpq(coeff, polyQ, targetIndex);
		
		// 4. Multiply by (n-1)! / (k-1)!		
		fmpz* numPtr = fmpq_numref(coeff);
		fmpz* denPtr = fmpq_denref(coeff);
		
		fmpz_mul(coeffNum, nFact, numPtr);
		fmpz_mul(coeffDen, kFact, denPtr);
		fmpz_tdiv_q(countOut, coeffNum, coeffDen);	
		
		fmpq_clear(coeff);
		fmpz_clear(coeffNum);
		fmpz_clear(coeffDen);
	}
		
	//Clear
	fmpq_poly_clear(poly);
	fmpq_poly_clear(polyQ);
	
}
*/
void stirling2_max_initial_lt(int n, int k, int m, int r, fmpz_t nFact, fmpz_t kFact, fmpq_mat_t coeffs, fmpz_t countOut) { //Q-nomial EGF
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
	fmpq_poly_t poly;
	fmpq_poly_t polyQ;
	fmpq_t prevCoeff;
	
	fmpq_poly_init(poly); 
	fmpq_poly_init(polyQ); 
	fmpq_init(prevCoeff);	
	fmpq_one(prevCoeff);

	for (int i = 0; i < m; i++) {
		fmpq* curCoeff = fmpq_mat_entry(coeffs, 0, i);
		fmpq_poly_set_coeff_fmpq(poly, i, curCoeff);		
		
		// Check if this power is valid
		if (i >= (r-1)) {			
			fmpq_poly_set_coeff_fmpq(polyQ, i, prevCoeff);
		}

		fmpq_set(prevCoeff, curCoeff);
	}
	
	fmpq_clear(prevCoeff);
	
	// Compute powers
	int targetIndex = n - k;
	if (k == 1) {
		fmpq_poly_one(poly);	
	}
	else {
		fmpq_poly_pow_trunc(poly, poly, k-1, targetIndex+1); 
	}
	//Multiply polynomials
	fmpq_poly_mullow(polyQ, polyQ, poly, targetIndex+1);
	
	//Extract coefficient
    if (targetIndex >= 0 && targetIndex < fmpq_poly_length(polyQ)) {
		fmpq_t coeff;
		fmpz_t coeffNum;
		fmpz_t coeffDen;
		
		fmpq_init(coeff);
		fmpz_init(coeffNum);
		fmpz_init(coeffDen);
		
		fmpq_poly_get_coeff_fmpq(coeff, polyQ, targetIndex);
		
		// 4. Multiply by (n-1)! / (k-1)!		
		fmpz* numPtr = fmpq_numref(coeff);
		fmpz* denPtr = fmpq_denref(coeff);
		
		fmpz_mul(coeffNum, nFact, numPtr);
		fmpz_mul(coeffDen, kFact, denPtr);
		fmpz_tdiv_q(countOut, coeffNum, coeffDen);	
		
		fmpq_clear(coeff);
		fmpz_clear(coeffNum);
		fmpz_clear(coeffDen);
	}
	
	
	//Clear
	fmpq_poly_clear(poly);
	fmpq_poly_clear(polyQ);
	
}
void stirling2_max_initial_ge(int n, int k, int m, int r, fmpz_t nFact, fmpz_t kFact, fmpq_mat_t coeffs, fmpz_t countOut) {
	//This calculates the Stirling2 number of set-parts with n elements, k parts, 
	//and having a max part size that is exactly m, and that has an initial part size greater-than-or-equal to r
	fmpz_t mCount;
	fmpz_init(mCount);
	stirling2_max_initial_lt(n, k, m, r, nFact, kFact, coeffs, countOut);
	stirling2_max_initial_lt(n, k, m-1, r, nFact, kFact, coeffs, mCount);	
	fmpz_sub(countOut, countOut, mCount);
	fmpz_clear(mCount);
}
void stirling2_max_initial_gt(int n, int k, int m, int r, fmpz_t nFact, fmpz_t kFact, fmpq_mat_t coeffs, fmpz_t countOut) {
	//This calculates the Stirling2 number of set-parts with n elements, k parts, 
	//and having a max part size that is exactly m, and that has an initial part size greater-than r
	fmpz_t mCount;
	fmpz_init(mCount);
	stirling2_max_initial_lt(n, k, m, r+1, nFact, kFact, coeffs, countOut);
	stirling2_max_initial_lt(n, k, m-1, r+1, nFact, kFact, coeffs, mCount);	
	fmpz_sub(countOut, countOut, mCount);
	fmpz_clear(mCount);
}



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


/* //OLD - NON-PRECACHED 
void stirling2_max_initial_lt_old(int n, int k, int m, int r, fmpz_t nFact, fmpz_t kFact, fmpz_t countOut) { //Q-nomial EGF
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
		fmpq_poly_pow_trunc(poly, poly, k-1, targetIndex+1); 
	}
	//Multiply polynomials
	fmpq_poly_mullow(polyQ, polyQ, poly, targetIndex+1);
	
	//Extract coefficient
    if (targetIndex >= 0 && targetIndex < fmpq_poly_length(polyQ)) {
		fmpq_poly_get_coeff_fmpq(coeff, polyQ, targetIndex);
		
		// 4. Multiply by (n-1)! / (k-1)!		
		fmpz* numPtr = fmpq_numref(coeff);
		fmpz* denPtr = fmpq_denref(coeff);
		
		fmpz_mul(coeffNum, nFact, numPtr);
		fmpz_mul(coeffDen, kFact, denPtr);
		fmpz_tdiv_q(countOut, coeffNum, coeffDen);	
	}
	
	
	//Clear
	fmpq_clear(coeff);
	fmpz_clear(coeffNum);
	fmpz_clear(coeffDen);
	fmpz_clear(prevDen);
	fmpq_poly_clear(poly);
	fmpq_poly_clear(polyQ);
	
}

void stirling2_max_lt_old(int n, int k, int m, fmpz_t nFact, fmpz_t kFact, fmpz_t countOut) { //Polynomial EGF
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
		fmpz* numPtr = fmpq_numref(coeff);
		fmpz* denPtr = fmpq_denref(coeff);
		
		fmpz_mul(coeffNum, nFact, numPtr);		
		fmpz_mul(coeffDen, kFact, denPtr);
		fmpz_tdiv_q(countOut, coeffNum, coeffDen);		
	}
	
	//Clear
	fmpq_clear(coeff);
	fmpz_clear(coeffNum);
	fmpz_clear(coeffDen);
	fmpq_poly_clear(poly);
}
*/