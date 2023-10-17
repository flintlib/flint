/*
    Copyright (C) 2007 William Hart and David Harvey
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(fmpz_poly_CRT_ui_unsigned, state)
{
    int i;

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

        fmpz_poly_randtest_unsigned(A, state, length, bits);

        fmpz_init(mod);
        num_primes = 0;
        primes[0] = n_nextprime(UWORD(1) << prime_bits, 0);
        fmpz_set_ui(mod, primes[0]);

        while (fmpz_bits(mod) <= bits)
        {
            primes[num_primes + 1] = n_nextprime(primes[num_primes], 0);
            fmpz_mul_ui(mod, mod, primes[num_primes + 1]);
            num_primes++;
        }

        num_primes++;

        nmod_poly_init(Amod, primes[0]);
        fmpz_poly_get_nmod_poly(Amod, A);
        fmpz_poly_set_nmod_poly_unsigned(B, Amod);
        fmpz_set_ui(mod, primes[0]);

        for (j = 1; j < num_primes; j++)
        {
            nmod_poly_clear(Amod);
            nmod_poly_init(Amod, primes[j]);
            fmpz_poly_get_nmod_poly(Amod, A);
            fmpz_poly_CRT_ui(B, B, mod, Amod, 0);
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
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(Amod);
        fmpz_poly_clear(A);
        fmpz_poly_clear(B);
        fmpz_poly_clear(C);
        fmpz_clear(mod);
    }

    TEST_FUNCTION_END(state);
}
