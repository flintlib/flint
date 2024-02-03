/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_factor_power235, state)
{
    slong ix;
    int result;

    for (ix = 0; ix < 3000 * flint_test_multiplier(); ix++)
    {
        ulong factor, exp, n1, n2, n1pow;
        int bits, type;

        type = n_randint(state, 4);

        if (type == 0)
        {
            /* Test random squares */
            bits = n_randint(state, FLINT_BITS / 2) + 1;
            n1 = n_randtest_bits(state, bits);
            n1pow = n1 * n1;
        }
        else if (type == 1)
        {
            /* Test random cubes */
            bits = n_randint(state, FLINT_BITS / 3) + 1;
            n1 = n_randtest_bits(state, bits);
            n1pow = n1 * n1 * n1;
        }
        else if (type == 2)
        {
            /* Test random fifth powers */
            bits = n_randint(state, FLINT_BITS / 5) + 1;
            n1 = n_randtest_bits(state, bits);
            n1pow = n1 * n1 * n1 * n1 * n1;
        }
        else
        {
            /* Test non 235-powers */
            do
                n1 = n_randtest(state);
            while (n_is_perfect_power235(n1));

            result = (!n_factor_power235(&exp, n1));
        }

        if (type < 3)
        {
            factor = n_factor_power235(&exp, n1pow);
            n2 = n_pow(factor, exp);
            result = (n1pow == n2);
        }

        if (!result)
            TEST_FUNCTION_FAIL(
                    "type %d\n"
                    "n1 = %wu, exp = %wu\n",
                    type, n1, exp);
    }

    TEST_FUNCTION_END(state);
}
