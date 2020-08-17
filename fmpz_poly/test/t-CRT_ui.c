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
#include "fmpz_poly.h"
#include "nmod_poly.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);

    flint_printf("CRT_ui....");
    fflush(stdout);

    

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong bits, prime_bits, length, num_primes, j;
        fmpz_t mod;
        fmpz_poly_t A, B, C;
        nmod_poly_t Amod;
        mp_limb_t primes[1000];

        bits = n_randint(state, 500) + 1;
        length = n_randint(state, 30) + 1;
        prime_bits = 1 + n_randint(state, FLINT_BITS - 1);

        fmpz_poly_init(A);
        fmpz_poly_init(B);
        fmpz_poly_init(C);

        fmpz_poly_randtest(A, state, length, bits);

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

        nmod_poly_init(Amod, primes[0]);
        fmpz_poly_get_nmod_poly(Amod, A);
        fmpz_poly_set_nmod_poly(B, Amod);
        fmpz_set_ui(mod, primes[0]);

        for (j = 1; j < num_primes; j++)
        {
            nmod_poly_clear(Amod);
            nmod_poly_init(Amod, primes[j]);
            fmpz_poly_get_nmod_poly(Amod, A);
            fmpz_poly_CRT_ui(B, B, mod, Amod, 1);
            fmpz_mul_ui(mod, mod, primes[j]);
        }

        if (!fmpz_poly_equal(B, A))
        {
            flint_printf("FAIL!\n");
            flint_printf("primes: ");
            for (j = 0; j < num_primes; j++)
                flint_printf("%wu ", primes[j]);
            flint_printf("\nA: \n");
            fmpz_poly_print(A);
            flint_printf("\nB: \n");
            fmpz_poly_print(B);
            flint_printf("\n");
            abort();
        }

        nmod_poly_clear(Amod);
        fmpz_poly_clear(A);
        fmpz_poly_clear(B);
        fmpz_poly_clear(C);
        fmpz_clear(mod);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
