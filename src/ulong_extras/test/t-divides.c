/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_divides, state)
{
    int i, result;

    for (i = 0; i < 2000 * flint_test_multiplier(); i++)
    {
        ulong n, p, q;
        int nbits, pbits;
        int flag, type;

        type = n_randint(state, 2);
        pbits = n_randint(state, FLINT_BITS + 1);

        if (type == 0)
        {
            /* test random values */
            nbits = n_randint(state, FLINT_BITS + 1);

            n = n_randtest_bits(state, nbits);
            p = n_randtest_bits(state, pbits);

            flag = n_divides(&q, n, p);

            result = ((flag && ((p == 0 && n == 0) || p*q == n)) ||
                    (!flag && q == 0 && ((p == 0 && n != 0) || p*q != n)));
        }
        else
        {
            /* test known divisible values */
            ulong s;
            int sbits = n_randint(state, FLINT_BITS - pbits + 1);

            p = n_randtest_bits(state, pbits);
            s = n_randtest_bits(state, sbits);

            n = p * s;

            flag = n_divides(&q, n, p);
            result = (flag && ((p == 0 && n == 0) || p*q == n));
        }

        if (!result)
            TEST_FUNCTION_FAIL("n = %wu, p = %wu, q = %wu\n", n, p, q);
    }

    TEST_FUNCTION_END(state);
}
