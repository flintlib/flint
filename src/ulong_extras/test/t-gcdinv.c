/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_gcdinv, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        ulong a, b, c, g, g2, s, t2, t, bits1, bits2, bits3, ainv;

        bits1 = n_randint(state, FLINT_BITS - 1) + 2;
        bits2 = n_randint(state, bits1) + 1;
        bits3 =
            bits1 == FLINT_BITS ? 0 : n_randint(state, FLINT_BITS - bits1) + 1;

        do
        {
            a = n_randtest_bits(state, bits1);
            b = n_randtest_bits(state, bits2);
        } while (n_gcd(a, b) != UWORD(1) || b >= a);

        c = bits3 == 0 ? 1 : n_randtest_bits(state, bits3);

        /* compare n_gcdinv with n_xgcd */
        g = n_xgcd(&s, &t, a * c, b * c);
        g2 = n_gcdinv(&t2, b * c, a * c);

        /* compute second cofactor modulo ac */
        t %= (a * c);           /* t is non-negative... */
        t = a * c - t;          /* ... but minus the actual cofactor */

        result = (g == g2 && t == t2);
        if (!result)
        {
            flint_printf("FAIL\n");
            flint_printf("Cofactor doesn't agree with n_xgcd\n");
            flint_printf("a = %wu, b = %wu, c = %wu\n", a, b, c);
            flint_printf("g = %wu, g2 = %wu, t = %wd, t2 = %wd\n", g, g2, t,
                         t2);
            fflush(stdout);
            flint_abort();
        }

        /* test b*t2 == 1 mod a */
        ainv = n_preinvert_limb(a);

        s = n_mulmod2_preinv(t2, b, a, ainv);

        result = (s == 1);
        if (!result)
        {
            flint_printf("FAIL\n");
            flint_printf("Incorrect inverse\n");
            flint_printf("a = %wu, b = %wu, c = %wu\n", a, b, c);
            flint_printf("g2 = %wu, s = %wd, t2 = %wd\n", g2, s, t2);
            fflush(stdout);
            flint_abort();
        }
    }

    /* test modulo 1 */
    {
        ulong s, g;

        g = n_gcdinv(&s, 0, 1);

        result = (g == 1 && s == 0);
        if (!result)
        {
            flint_printf("FAIL\n");
            flint_printf("Incorrect modulo 1\n");
            flint_printf("g = %wu, s = %wu\n", g, s);
            fflush(stdout);
            flint_abort();
        }
    }

    /* check gcd not 1 when a = 0 (and b != 1) */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        ulong b, s, g;

        b = n_randtest_not_zero(state);

        g = n_gcdinv(&s, 0, b);

        result = (g != 1 || b == 1);
        if (!result)
        {
            flint_printf("FAIL\n");
            flint_printf("gcd(0, b) == 1\n");
            flint_printf("b = %wu, s = %wu\n", b, s);
            fflush(stdout);
            flint_abort();
        }
    }

    TEST_FUNCTION_END(state);
}
