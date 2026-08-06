/*
    Copyright (C) 2009, 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_divrem2_preinv, state)
{
    int i, result;

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        ulong d, dinv, n, q1, q2, r1, r2, norm;

        d = n_randtest_not_zero(state);
        n = n_randtest(state);

        norm = flint_clz(d);
        dinv = n_preinvert_limb_prenorm(d << norm);

        switch (i % 3)
        {
            case 0:
                r1 = n_divrem2_preinv(&q1, n, d, dinv);
                break;
            case 1:
                if (norm == 0)
                    r1 = n_divrem_norm(&q1, n, d);
                else
                    r1 = n_divrem_preinv_unnorm(&q1, n, d, dinv, norm);
                break;
            default:
                r1 = n_divrem_preinv(&q1, n, d, dinv, norm);
                break;
        }

        q2 = n / d;
        r2 = n % d;

        result = (q1 == q2 && r1 == r2);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "i = %d, n = %wu, d = %wu, dinv = %wu\n"
                    "q1 = %wu, q2 = %wu\n"
                    "r1 = %wu, r2 = %wu\n",
                    i, n, d, dinv, q1, q2, r1, r2);
    }

    TEST_FUNCTION_END(state);
}
