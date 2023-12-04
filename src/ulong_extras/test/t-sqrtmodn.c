/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_sqrtmodn, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test random squares mod n */
    {
        mp_limb_t a, b, n, ninv;
        slong num, i;
        flint_bitcnt_t bits;
        mp_limb_t * sqrt;
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
        {
            flint_printf("FAIL:\n");
            flint_printf("n = %wu\n", n);
            flint_printf("a = %wu\n", a);
            flint_printf("b = %wu\n", b);
            flint_printf("num = %wd\n", num);

            if (!btest)
                flint_printf("Square root not found.\n");
            if (i != num)
                flint_printf("%wu not a square root of %wu mod %wu\n", sqrt[i], a, n);

            fflush(stdout);
            flint_abort();
        }

        flint_free(sqrt);
    }

    for (i = 0; i < 500 * flint_test_multiplier(); i++) /* test random nonsquares */
    {
        mp_limb_t a, b, n, ninv;
        flint_bitcnt_t bits;
        mp_limb_t * sqrt;
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
            {
                flint_printf("FAIL:\n");
                flint_printf("%wu^2 is not %wu mod %wu\n", sqrt[0], a, n);
                fflush(stdout);
                flint_abort();
            }

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
        {
            flint_printf("FAIL:\n");
            flint_printf("n = %wu\n", n);
            flint_printf("a = %wu\n", a);
            flint_printf("b = %wu\n", b);

            fflush(stdout);
            flint_abort();
        }

        flint_free(sqrt);
    }

    TEST_FUNCTION_END(state);
}
