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

TEST_FUNCTION_START(n_mulmod2, state)
{
    int i, result;

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        ulong a, b, d, r1, r2, q, p1, p2;

        d = n_randtest_not_zero(state);
        a = n_randtest(state) % d;
        b = n_randtest(state) % d;

        r1 = n_mulmod2(a, b, d);

        umul_ppmm(p1, p2, a, b);
        p1 %= d;
        udiv_qrnnd(q, r2, p1, p2, d);

        result = (r1 == r2);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = %wu, b = %wu, d = %wu\n", a, b, d);
            flint_printf("q = %wu, r1 = %wu, r2 = %wu\n", q, r1, r2);
            fflush(stdout);
            flint_abort();
        }
    }

    TEST_FUNCTION_END(state);
}
