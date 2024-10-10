/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_sqrtmodn, state)
{
    slong ix;
    int result;

    for (ix = 0; ix < 1000 * flint_test_multiplier(); ix++) /* Test random squares mod n */
    {
        ulong a, b, n, ninv;
        slong num, i;
        flint_bitcnt_t bits;
        ulong * sqrt;
        int btest;
        n_factor_t fac;

        bits = n_randint(state, 18) + 1;
        n = n_randtest_bits(state, bits);
        b = n_randtest(state) % n;

        n_factor_init(&fac);
        n_factor(&fac, n, 0);

        ninv = n_preinvert_limb(n);
        a = n_mulmod2_preinv(b, b, n, ninv);

        num = n_sqrtmodn(&sqrt, a, &fac);

        btest = 0;
        for (i = 0; i < num; i++)
        {
            if (a != n_mulmod2_preinv(sqrt[i], sqrt[i], n, ninv))
                break;
            if (sqrt[i] == b)
                btest = 1;
        }

        result = btest & (i == num);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "n = %wu\n"
                    "a = %wu\n"
                    "b = %wu\n"
                    "num = %wd\n"
                    "btest = %d\n"
                    "i != num = %d\n",
                    n, a, b, num, btest, i != num);

        flint_free(sqrt);
    }

    for (ix = 0; ix < 500 * flint_test_multiplier(); ix++) /* test random nonsquares */
    {
        ulong a, b, n, ninv;
        flint_bitcnt_t bits;
        ulong * sqrt;
        n_factor_t fac;

        bits = n_randint(state, 18) + 2;
        n = n_randtest_bits(state, bits);
        if (n == 2) n++;
        n_factor_init(&fac);
        n_factor(&fac, n, 0);

        ninv = n_preinvert_limb(n);

        a = n_randtest(state) % n;
        while (n_sqrtmodn(&sqrt, a, &fac))
        {
            if (n_mulmod2_preinv(sqrt[0], sqrt[0], n, ninv) != a)
                TEST_FUNCTION_FAIL("%wu^2 is not %wu mod %wu\n", sqrt[0], a, n);

            flint_free(sqrt);
            a = n_randtest(state) % n;
        }

        for (b = 0; b < n; b++)
        {
            if (n_mulmod2_preinv(b, b, n, ninv) == a)
                break;
        }

        result = (b == n);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "n = %wu\n"
                    "a = %wu\n"
                    "b = %wu\n",
                    n, a, b);

        flint_free(sqrt);
    }

    TEST_FUNCTION_END(state);
}
