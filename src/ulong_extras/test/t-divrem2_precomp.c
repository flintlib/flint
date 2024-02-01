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

TEST_FUNCTION_START(n_divrem2_precomp, state)
{
    int result;
    ulong i;

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        mp_limb_t d, n, r1, r2, q1, q2;
        double dpre;

        d = n_randtest_not_zero(state);
        n = n_randtest(state);

        dpre = n_precompute_inverse(d);

        r1 = n_divrem2_precomp(&q1, n, d, dpre);
        r2 = n % d;
        q2 = n / d;

        result = ((r1 == r2) && (q1 == q2));
        if (!result)
            TEST_FUNCTION_FAIL(
                    "n = %wu, d = %wu, dpre = %f\n"
                    "q1 = %wu, q2 = %wu, r1 = %wu, r2 = %wu\n",
                    n, d, dpre, q1, q2, r1, r2);
    }

    TEST_FUNCTION_END(state);
}
