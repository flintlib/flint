/*
    Copyright (C) 2009, 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_gcd, state)
{
    int i, result;

    /* test gcd(ac, bc) == gcd(a, b) */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        ulong a, b, c, g, bits1, bits2, bits3, mbits;

        bits1 = n_randint(state, FLINT_BITS - 1) + 1;
        bits2 = n_randint(state, FLINT_BITS - 1) + 1;
        mbits = FLINT_MAX(bits1, bits2);

        bits3 = mbits == FLINT_BITS ?
            0 : n_randint(state, FLINT_BITS - mbits) + 1;

        do
        {
            a = n_randtest_bits(state, bits1);
            b = n_randtest_bits(state, bits2);
        } while (n_gcd(a, b) != UWORD(1));

        c = bits3 == 0 ? 1 : n_randtest_bits(state, bits3);

        g = n_gcd(a * c, b * c);

        result = (g == c);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("gcd(ac, bc) != gcd(a, b)\n");
            flint_printf("a = %wu, b = %wu, c = %wu, g = %wu\n", a, b, c, g);
            fflush(stdout);
            flint_abort();
        }
    }

    /* test gcd(a, 0) == a */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        ulong a, g;

        a = n_randtest(state);

        g = n_gcd(a, 0);

        result = (g == a);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("gcd(a, 0) != a\n");
            flint_printf("a = %wu\n", a);
            fflush(stdout);
            flint_abort();
        }
    }

    /* test gcd(0, b) == b */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        ulong b, g;

        b = n_randtest(state);

        g = n_gcd(0, b);

        result = (g == b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("gcd(0, b) != b\n");
            flint_printf("b = %wu\n", b);
            fflush(stdout);
            flint_abort();
        }
    }

    /* test gcd(0, 0) == 0 */
    {
        result = (n_gcd(0, 0) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("gcd(0, 0) != 0\n");
            fflush(stdout);
            flint_abort();
        }
    }

    TEST_FUNCTION_END(state);
}
