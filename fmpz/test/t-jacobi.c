/*
    Copyright (C) 2009 William Hart

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
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("jacobi....");
    fflush(stdout);

    
    _flint_rand_init_gmp(state);

    for (i = 0; i < 3000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, p;
        mpz_t b, q;
        int r1, r2;

        fmpz_init(a);
        fmpz_init(p);

        mpz_init(b);
        mpz_init(q);

        mpz_rrandomb(q, state->gmp_state, n_randint(state, 200) + 1);
#ifdef mpz_next_likely_prime
        mpz_next_likely_prime(q, q, state->gmp_state);
#else
        mpz_nextprime(q, q);
#endif
        fmpz_set_mpz(p, q);

        mpz_rrandomb(b, state->gmp_state, n_randint(state, 200) + 1);
        mpz_mod(b, b, q);
        if (n_randint(state, 2))
            mpz_neg(b, b);
        fmpz_set_mpz(a, b);

        r1 = fmpz_jacobi(a, p);
        r2 = mpz_jacobi(b, q);
        result = (r1 == r2);

        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("b = %Zd, q = %Zd\n", b, q);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(p);

        mpz_clear(b);
        mpz_clear(q);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
