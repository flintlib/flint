/*
    Copyright (C) 2009, 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_xgcd, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        ulong a, b, c, g, bits1, bits2, bits3, ph, pl, qh, ql;
        ulong s, t;

        bits1 = n_randint(state, FLINT_BITS - 1) + 2;
        bits2 = n_randint(state, bits1) + 1;
        bits3 =
            bits1 == FLINT_BITS ? 0 : n_randint(state, FLINT_BITS - bits1) + 1;

        do
        {
            a = n_randtest_bits(state, bits1);
            b = n_randtest_bits(state, bits2);
        } while (n_gcd(a, b) != UWORD(1) || b > a);

        c = bits3 == 0 ? 1 : n_randtest_bits(state, bits3);

        g = n_xgcd(&s, &t, a * c, b * c);

        /* s*ac + t*bc */
        umul_ppmm(ph, pl, a * c, s);
        umul_ppmm(qh, ql, b * c, t);
        sub_ddmmss(ph, pl, ph, pl, qh, ql);

        result = (s <= b * c || b * c == 0);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "s >= b*c\n"
                    "a = %wu, b = %wu, c = %wu, g = %wu, s = %wu, t = %wu\n",
                    a, b, c, g, s, t);

        result = (t <= a * c || a * c == 0);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "t >= a*c\n"
                    "a = %wu, b = %wu, c = %wu, g = %wu, s = %wu, t = %wu\n",
                    a, b, c, g, s, t);

        result = (g == c && ph == UWORD(0) && pl == c);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "g != c or s*ac + t*bc != c\n"
                    "a = %wu, b = %wu, c = %wu, g = %wu, s = %wu, t = %wu\n",
                    a, b, c, g, s, t);
    }

    /* a = 0, b = 0 */
    {
        ulong g, s, t;

        g = n_xgcd(&s, &t, 0, 0);

        result = (g == 0 && s == 1 && t == 0);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "Case a = 0, b = 0\n"
                    "g = %wu, s = %wu, t = %wu\n",
                    g, s, t);
    }

    /* b = 0 */
    {
        ulong a, g, s, t;

        a = n_randtest(state);

        g = n_xgcd(&s, &t, a, 0);

        result = (g == a && s == 1 && t == 0);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "Case a = 0\n"
                    "a = %wu, g = %wu, s = %wu, t = %wu\n",
                    a, g, s, t);
    }

    TEST_FUNCTION_END(state);
}
