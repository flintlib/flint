#include <stdio.h>
#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

int main()
{
    long i;
    fmpz_t x, y, prod;

    /* Data needed by multi CRT functions */
    fmpz_comb_t comb;
    fmpz ** comb_temp;
    fmpz_t temp_fmpz, temp_fmpz2;
    mp_limb_t * primes;
    mp_limb_t * residues;

    long num_primes;

    fmpz_init(x);
    fmpz_init(y);
    fmpz_init(prod);

    fmpz_set_str(x, "-12345678901234567890", 10);

    num_primes = 18;   /* We assume the number of needed primes is known */

    primes = malloc(num_primes * sizeof(mp_limb_t));
    residues = malloc(num_primes * sizeof(mp_limb_t));

    primes[0] = 2;
    for (i = 1; i < num_primes; i++)
        primes[i] = n_nextprime(primes[i-1], 0);

    fmpz_comb_init(comb, primes, num_primes);
    comb_temp = fmpz_comb_temp_init(comb);
    fmpz_init(temp_fmpz);
    fmpz_init(temp_fmpz2);

    /* Reduce modulo all primes */
    fmpz_multi_mod_ui(residues, x, comb, comb_temp, temp_fmpz);

    /* Reconstruct */
    fmpz_multi_CRT_ui(y, residues, comb, comb_temp, temp_fmpz, temp_fmpz2);

    for (i = 0; i < num_primes; i++)
        printf("residue mod %lu = %lu\n", residues[i], primes[i]);

    printf("reconstruction = ");
    fmpz_print(y);
    printf("\n");

    fmpz_clear(x);
    fmpz_clear(y);

    fmpz_comb_temp_free(comb, comb_temp);
    fmpz_comb_clear(comb);
    fmpz_clear(temp_fmpz);
    fmpz_clear(temp_fmpz2);

    free(residues);
    free(primes);
}
