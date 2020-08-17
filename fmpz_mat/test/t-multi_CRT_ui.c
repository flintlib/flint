/*
    Copyright (C) 2007 William Hart and David Harvey
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "nmod_mat.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);

    flint_printf("multi_CRT_ui....");
    fflush(stdout);

    

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong bits, prime_bits, rows, cols, num_primes, j;
        fmpz_t mod;
        fmpz_mat_t A, B, C;
        nmod_mat_t Amod[1000];
        mp_limb_t primes[1000];

        bits = n_randint(state, 500) + 1;
        rows = n_randint(state, 10);
        cols = n_randint(state, 10);
        prime_bits = 1 + n_randint(state, FLINT_BITS - 1);

        fmpz_mat_init(A, rows, cols);
        fmpz_mat_init(B, rows, cols);
        fmpz_mat_init(C, rows, cols);

        fmpz_mat_randtest(A, state, bits);

        fmpz_init(mod);
        num_primes = 0;
        primes[0] = n_nextprime(UWORD(1) << prime_bits, 0);
        fmpz_set_ui(mod, primes[0]);

        /* + 1 for sign */
        while (fmpz_bits(mod) <= bits + 1)
        {
            primes[num_primes + 1] = n_nextprime(primes[num_primes], 0);
            fmpz_mul_ui(mod, mod, primes[num_primes + 1]);
            num_primes++;
        }

        num_primes++;

        for (j = 0; j < num_primes; j++)
            nmod_mat_init(Amod[j], rows, cols, primes[j]);

        fmpz_mat_multi_mod_ui(Amod, num_primes, A);
        fmpz_mat_multi_CRT_ui(B, Amod, num_primes, 1);

        if (!fmpz_mat_equal(B, A))
        {
            flint_printf("FAIL!\n");
            flint_printf("primes: ");
            for (j = 0; j < num_primes; j++)
                flint_printf("%wu ", primes[j]);
            flint_printf("\nA: \n");
            fmpz_mat_print_pretty(A);
            flint_printf("\nB: \n");
            fmpz_mat_print_pretty(B);
            flint_printf("\n");
            abort();
        }

        for (j = 0; j < num_primes; j++)
            nmod_mat_clear(Amod[j]);
        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_clear(mod);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
