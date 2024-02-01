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

TEST_FUNCTION_START(n_factor, state)
{
    int i, j, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        mp_limb_t n1, n2;
        n_factor_t factors;
        int type;

        n_factor_init(&factors);

        type = n_randint(state, 2);

        if (type == 0)
        {
            /* Test random numbers */
            n1 = n_randtest_not_zero(state);
        }
        else
        {
            /* Test random prime */
            n1 = n_randprime(state, 2 + n_randint(state, FLINT_BITS - 1), 1);
        }

        n_factor(&factors, n1, 0);

        n2 = UWORD(1);
        for (j = 0; j < factors.num; j++)
        {
            n2 *= n_pow(factors.p[j], factors.exp[j]);
        }

        result = (n1 == n2);
        if (!result)
            TEST_FUNCTION_FAIL("n1 = %wu, n2 = %wu\n", n1, n2);
    }

    TEST_FUNCTION_END(state);
}
