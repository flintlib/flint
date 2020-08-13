/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Demo FLINT program for balanced multimodular reduction and
    reconstruction using the Chinese Remainder Theorem.
*/

#include <stdlib.h>
#include <stdio.h>
#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

int main(int argc, char* argv[])
{
    slong i;
    fmpz_t x, y;

    /* Data needed by multi CRT functions */
    fmpz_comb_t comb;
    fmpz_comb_temp_t comb_temp;

    mp_limb_t * primes;
    mp_limb_t * residues;

    slong num_primes;

    if (argc != 3)
    {
        flint_printf("Syntax: crt <integer> <num_primes>\n");
        return EXIT_FAILURE;
    }

    num_primes = atoi(argv[2]);

    if (num_primes < 1)
    {
        flint_printf("Requires num_primes >= 1\n");
        return EXIT_FAILURE;
    }

    fmpz_init(x);
    fmpz_init(y);

    fmpz_set_str(x, argv[1], 10);

    primes = flint_malloc(num_primes * sizeof(mp_limb_t));
    residues = flint_malloc(num_primes * sizeof(mp_limb_t));

    primes[0] = 2;
    for (i = 1; i < num_primes; i++)
        primes[i] = n_nextprime(primes[i-1], 0);

    fmpz_comb_init(comb, primes, num_primes);
    fmpz_comb_temp_init(comb_temp, comb);

    /* Reduce modulo all primes */
    fmpz_multi_mod_ui(residues, x, comb, comb_temp);

    /* Reconstruct */
    fmpz_multi_CRT_ui(y, residues, comb, comb_temp, 1);

    for (i = 0; i < num_primes; i++)
        flint_printf("residue mod %wu = %wu\n", primes[i], residues[i]);

    flint_printf("reconstruction = ");
    fmpz_print(y);
    flint_printf("\n");

    fmpz_clear(x);
    fmpz_clear(y);

    fmpz_comb_temp_clear(comb_temp);
    fmpz_comb_clear(comb);

    flint_free(residues);
    flint_free(primes);

    return EXIT_SUCCESS;
}
