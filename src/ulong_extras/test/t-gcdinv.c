/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_gcdinv, state)
{
    slong ix;
    int result;

    /* test modulo 1 */
    {
        ulong a, b, s, g;

        a = 0;
        b = 1;
        g = n_gcdinv(&s, a, b);

        result = (g == 1 && s == 0);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "GCD modulo 1 does not return g = 1 and s = 0\n"
                    "g = %wu, s = %wu\n",
                    g, s);
    }

    /* check gcd not 1 when a = 0 (and b != 1) */
    for (ix = 0; ix < 1000 * flint_test_multiplier(); ix++)
    {
        ulong a, b, s, g;

        a = 0;
        b = n_randtest_not_zero(state);
        g = n_gcdinv(&s, a, b);

        result = (g != 1 || b == 1);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "gcd(0, b) == 1\n"
                    "b = %wu, s = %wu\n",
                    b, s);
    }

    for (ix = 0; ix < 10000 * flint_test_multiplier(); ix++)
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
            TEST_FUNCTION_FAIL(
                    "Cofactor doesn't agree with n_xgcd\n"
                    "a = %wu, b = %wu, c = %wu\n"
                    "g = %wu, g2 = %wu, t = %wd, t2 = %wd\n",
                    a, b, c, g, g2, t, t2);

        /* test b*t2 == 1 mod a */
        ainv = n_preinvert_limb(a);

        s = n_mulmod2_preinv(t2, b, a, ainv);

        result = (s == 1);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "Incorrect inverse\n"
                    "a = %wu, b = %wu, c = %wu\n"
                    "g2 = %wu, s = %wd, t2 = %wd\n",
                    a, b, c, g2, s, t2);
    }

    TEST_FUNCTION_END(state);
}
