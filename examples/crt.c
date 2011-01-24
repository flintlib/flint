#include <stdio.h>
#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

int main()
{
    long i, bit_bound;
    mp_limb_t prime, res;
    fmpz_t x, y, prod;

    fmpz_init(x);
    fmpz_init(y);
    fmpz_init(prod);

    fmpz_set_str(x, "-12345678901234567890", 10);
    bit_bound = fmpz_bits(x) + 2;

    fmpz_set_ui(y, 0);
    fmpz_set_ui(prod, 1);

    prime = 0;
    for (i = 0; fmpz_bits(prod) < bit_bound; i++)
    {
        prime = n_nextprime(prime, 0);
        res = fmpz_fdiv_ui(x, prime);
        fmpz_CRT_ui(y, y, prod, res, prime);

        printf("residue mod %lu = %lu; reconstruction = ", res, prime);
        fmpz_print(y);
        printf("\n");

        fmpz_mul_ui(prod, prod, prime);
    }

    fmpz_clear(x);
    fmpz_clear(y);
    fmpz_clear(prod);
}
