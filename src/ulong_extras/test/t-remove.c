/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_remove, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++) /* Test random numbers */
    {
        ulong n1, n2, orig_n;
        mpz_t d_n2, d_n1, d_p;
        int exp1, exp2;
        ulong j;

        mpz_init(d_n1);
        mpz_init(d_n2);
        mpz_init(d_p);

        n1 = n_randtest_not_zero(state);
        orig_n = n1;

        for (j = 0; j < FLINT_NUM_PRIMES_SMALL/10; j++)
        {
            flint_mpz_set_ui(d_n1, n1);
            flint_mpz_set_ui(d_p, flint_primes_small[j]);
            exp1 = n_remove(&n1, flint_primes_small[j]);
            exp2 = mpz_remove(d_n2, d_n1, d_p);
            n2 = flint_mpz_get_ui(d_n2);

            result = ((exp1 == exp2) && (n1 == n2));
            if (!result)
                TEST_FUNCTION_FAIL("n = %wu, exp1 = %d, exp2 = %d, n1 = %wu, n2 = %wu, p = %d\n", orig_n, exp1, exp2, n1, n2, flint_primes_small[j]);
        }

        mpz_clear(d_n1);
        mpz_clear(d_n2);
        mpz_clear(d_p);
    }

    for (i = 0; i < 10000 * flint_test_multiplier(); i++) /* Test perfect powers */
    {
        ulong n1, n2, orig_n, base;
        mpz_t d_n2, d_n1, d_p;
        int exp1, exp2, exp;
        ulong j;

        mpz_init(d_n1);
        mpz_init(d_n2);
        mpz_init(d_p);

        base = n_randtest_not_zero(state);
        n1 = base;
        exp = n_randint(state, FLINT_BITS/FLINT_BIT_COUNT(n1)) + 1;
        n1 = n_pow(base, exp);

        orig_n = n1;

        for (j = 0; j < FLINT_NUM_PRIMES_SMALL/10; j++)
        {
            flint_mpz_set_ui(d_n1, n1);
            flint_mpz_set_ui(d_p, flint_primes_small[j]);
            exp1 = n_remove(&n1, flint_primes_small[j]);
            exp2 = mpz_remove(d_n2, d_n1, d_p);
            n2 = flint_mpz_get_ui(d_n2);

            result = ((exp1 == exp2) && (n1 == n2));
            if (!result)
                TEST_FUNCTION_FAIL("n = %wu, exp1 = %d, exp2 = %d, n1 = %wu, n2 = %wu, p = %d\n", orig_n, exp1, exp2, n1, n2, flint_primes_small[j]);
        }

        mpz_clear(d_n1);
        mpz_clear(d_n2);
        mpz_clear(d_p);
    }

    for (i = 0; i < 10000 * flint_test_multiplier(); i++) /* Test arbitrary p, including p^2 >= 2^FLINT_BITS */
    {
        ulong n1, n2, orig_n, p, hi, lo;
        mpz_t d_n2, d_n1, d_p;
        int exp1, exp2;

        mpz_init(d_n1);
        mpz_init(d_n2);
        mpz_init(d_p);

        /* p^2 is 1, 1 and 0 modulo 2^FLINT_BITS */
        if (i == 0)
            p = UWORD_MAX;
        else if (i == 1)
            p = (UWORD(1) << (FLINT_BITS - 1)) + 1;
        else if (i == 2)
            p = UWORD(1) << (FLINT_BITS / 2);
        else
        {
            do
                p = n_randtest(state);
            while (p < 2);
        }

        /* n1 = r * p^k, multiplying in factors p until a random stop or overflow */
        n1 = (i < 3) ? p : n_randint(state, 1000) + 1;
        while (n_randint(state, 4) != 0)
        {
            umul_ppmm(hi, lo, n1, p);
            if (hi != 0)
                break;
            n1 = lo;
        }

        orig_n = n1;
        flint_mpz_set_ui(d_n1, n1);
        flint_mpz_set_ui(d_p, p);
        exp1 = n_remove(&n1, p);
        exp2 = mpz_remove(d_n2, d_n1, d_p);
        n2 = flint_mpz_get_ui(d_n2);

        result = ((exp1 == exp2) && (n1 == n2));
        if (!result)
            TEST_FUNCTION_FAIL("n = %wu, exp1 = %d, exp2 = %d, n1 = %wu, n2 = %wu, p = %wu\n", orig_n, exp1, exp2, n1, n2, p);

        mpz_clear(d_n1);
        mpz_clear(d_n2);
        mpz_clear(d_p);
    }

    TEST_FUNCTION_END(state);
}
