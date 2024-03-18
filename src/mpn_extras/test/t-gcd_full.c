/*
    Copyright (C) 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"

TEST_FUNCTION_START(flint_mpn_gcd_full, state)
{
    int i, result;
    mpz_t a, b, c, g;
    gmp_randstate_t st;
    slong s1, s2;

    mpz_init(a);
    mpz_init(b);
    mpz_init(c);
    /* don't init g */
    gmp_randinit_default(st);

    for (i = 0; i < 10000; i++)
    {
        do {
            mpz_urandomb(a, st, n_randint(state, 200));
        } while (mpz_sgn(a) == 0);
        do {
            mpz_urandomb(b, st, n_randint(state, 200));
        } while (mpz_sgn(b) == 0);
        do {
            mpz_urandomb(c, st, n_randint(state, 200));
        } while (mpz_sgn(c) == 0);

        mpz_mul(a, a, c);
        mpz_mul(b, b, c);
        mpz_mul_2exp(a, a, n_randint(state, 200));
        mpz_mul_2exp(b, b, n_randint(state, 200));

        mpz_gcd(c, a, b);

        s1 = (mpz_sizeinbase(a, 2) - 1)/FLINT_BITS + 1;
        s2 = (mpz_sizeinbase(b, 2) - 1)/FLINT_BITS + 1;

        g->_mp_d = flint_malloc(FLINT_MIN(s1, s2)*sizeof(mp_limb_t));

        g->_mp_size = flint_mpn_gcd_full(g->_mp_d, a->_mp_d, a->_mp_size, b->_mp_d, b->_mp_size);

        result = (mpz_cmp(g, c) == 0);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "%{mpz}\n"
                    "%{mpz}\n",
                    g, c);

        flint_free(g->_mp_d);
    }

    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(c);
    /* don't clear g */
    gmp_randclear(st);

    TEST_FUNCTION_END(state);
}
