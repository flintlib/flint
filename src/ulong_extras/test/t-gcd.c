/*
    Copyright (C) 2009, 2015 William Hart
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_gcd, state)
{
    int i, result;

    if (n_gcd(0, 0) != 0)
        TEST_FUNCTION_FAIL("gcd(0, 0) != 0\n");

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        ulong a, b, c, g;
        int type;

        type = n_randint(state, 30);

        if (type == 0)
        {
            /* gcd(a, 0) == a */
            a = n_randtest(state);
            b = 0;
            c = 0;
            g = n_gcd(a, b);

            result = (g == a);
        }
        else if (type == 1)
        {
            /* gcd(0, b) == b */
            a = 0;
            b = n_randtest(state);
            c = 0;
            g = n_gcd(a, b);

            result = (g == b);
        }
        else
        {
            /* gcd(ac, bc) == gcd(a, b) */
            int bits1, bits2, bits3, mbits;

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
        }

        if (!result)
            TEST_FUNCTION_FAIL(
                    "type %d\n"
                    "a = %wu, b = %wu, c = %wu, g = %wu\n",
                    type, a, b, c, g);
    }

    TEST_FUNCTION_END(state);
}
