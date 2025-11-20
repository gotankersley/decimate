#include <stdio.h>
#include <gmp.h>

//gcc gmp_test.c -o gmp_test -lgmp
int main() {
	mpz_t num;
	mpz_init_set_ui(num, 12345);
	gmp_printf("GMP number: %Zd\n", num);
	mpz_clear(num);
	return 0;
}