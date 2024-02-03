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

TEST_FUNCTION_START(n_factor_trial, state)
{
    slong ix, jx;
    int result;

    for (ix = 0; ix < 1000 * flint_test_multiplier(); ix++) /* Test random numbers */
    {
        ulong n1, n2;
        n_factor_t factors;

        n_factor_init(&factors);

        n1 = n_randtest_not_zero(state);
        n2 = n_factor_trial(&factors, n1, UWORD(10000));

        for (jx = 0; jx < factors.num; jx++)
        {
            n2 *= n_pow(factors.p[jx], factors.exp[jx]);
        }

        result = (n1 == n2);
        if (!result)
            TEST_FUNCTION_FAIL("n1 = %wu, n2 = %wu\n", n1, n2);
    }

    TEST_FUNCTION_END(state);
}
