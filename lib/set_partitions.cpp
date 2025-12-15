#include "set_partitions.h"
#include <iostream>

using std::cout, std::endl;

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

//Custom function to keep the numerators and denominators separate
void mult_poly(fmpz_mat_t coeffs1, fmpz_mat_t coeffs2, fmpz_mat_t coeffsOut) {
	//Based on:
	//def mult_poly(a, b):
	//	p = [0] * (len(a) + len(b) - 1)	
	//	for i, x in enumerate(a):
	//		for j, y in enumerate(b):			
	//			p[i + j] += x * y	
	//	return p
	int polyLen1 = coeffs1->c;
	int polyLen2 = coeffs2->c;	
	int resLen = (polyLen1 + polyLen2 - 1);
	fmpz_mat_init(coeffsOut, FRACTION_SIZE, resLen);
	fmpz_mat_zero(coeffsOut);
	for (int i = 0; i < resLen; i++) {
		fmpz_one(fmpz_mat_entry(coeffsOut, DEN, i));
	}
	
	fmpz_t productN;
	fmpz_t productD;
	
	fmpz_init(productN);
	fmpz_init(productD);
	
	for (int i = 0; i < polyLen1; i++) {
		fmpz* n1 = fmpz_mat_entry(coeffs1, NUM, i);
		fmpz* d1 = fmpz_mat_entry(coeffs1, DEN, i);
		
		for (int j = 0; j < polyLen2; j++) {			
			fmpz* n2 = fmpz_mat_entry(coeffs2, NUM, j);
			fmpz* d2 = fmpz_mat_entry(coeffs2, DEN, j);
		
			//First multiply the fractions
			fmpz_mul(productN, n1, n2);
			fmpz_mul(productD, d1, d2);
			
			//Second add the fractions, (have to make sure they have a common denominator)
			int idx = i+j;
			fmpz* nOut = fmpz_mat_entry(coeffsOut, NUM, idx);
			fmpz* dOut = fmpz_mat_entry(coeffsOut, DEN, idx);
			
			if (fmpz_equal(dOut, productD)) {
				//Already common denominator, can just add
				fmpz_add(nOut, nOut, productN); 
			}
			else {
				fmpz_mul(nOut, nOut, productD); 
				fmpz_mul(productN, productN, dOut); 
				fmpz_add(nOut, nOut, productN); 
				fmpz_mul(dOut, dOut, productD);
			}
		}
	}
	
	fmpz_clear(productN);
	fmpz_clear(productD);
	cout << "Poly mult res:" << endl;
	fmpz_mat_print_pretty(coeffsOut);
	cout << endl;
}

void gen_coeff_table(int maxN, int maxK, int m) {//Output is serialized
	//Note:  This is the precalculated part that is used by the
	//stirling2_max_less_than() function to make it viable.
	//We are calculating everything but the final step where n! is added
	fmpz_mat_t baseCoeffs;	
	gen_initial_coeffs(m, baseCoeffs);	
	//serialize_mat("K1.mat", baseCoeffs);
	cout << "Base Coeffs: " << endl;
	fmpz_mat_print_pretty(baseCoeffs);
	cout << endl;
	
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
		cout << "Col: " << k << endl;
		fmpz_mat_print_pretty(col);
		cout << endl;
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

void stirling2_max(fmpz_mat_t coeffs, int n, int k, int m, fmpz_t countOut) {
	//This calculates the Stirling2 number of set-parts with n elements, k parts, and having a max part size that is exactly m
	fmpz_t count2;
	fmpz_init(count2);
	stirling2_max_less_than_coeffs(coeffs, n, k, m-1, count2);	
	stirling2_max_less_than_coeffs(coeffs, n, k, m, countOut);
	fmpz_sub(countOut, countOut, count2);
	fmpz_clear(count2);
}

void stirling2_max_r(fmpz_mat_t coeffs, int n, int k, int m, int r, fmpz_t countOut) {
	//This calculates the Stirling2 number of set-parts with n elements, k parts, and having a max part size that is exactly m, and has an initial part size of r
	
	//TODO - See if this can be optimized
	if (r == m) { //This is the total minus all the other ones
		fmpz_t previousSum;
		fmpz_init(previousSum);
		fmpz_t count2;
		fmpz_init(count2);
		for (int i = 0; i < r; i++) {
			stirling2_max_r(coeffs, n, k, m, i, count2);
			fmpz_add(previousSum, previousSum, count2);
		}
		fmpz_clear(count2);
		stirling2_max_less_than_coeffs(coeffs, n, k, m, countOut); //Total
		fmpz_add(countOut, countOut, previousSum);
		fmpz_clear(previousSum);
	}
	else {		
		fmpz_t count2;
		fmpz_init(count2);
		stirling2_max_less_than_coeffs(coeffs, n-r, k-1, m-1, count2);
		stirling2_max_less_than_coeffs(coeffs, n-r, k-1, m, countOut);
		fmpz_sub(countOut, countOut, count2);
		fmpz_clear(count2);
		uint64_t cmb = comb(n-1, r-1);
		fmpz_mul_ui(countOut, countOut, cmb);		
	}
}


