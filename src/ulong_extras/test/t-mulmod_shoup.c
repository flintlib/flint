/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2015 Vladimir Glazachev
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_mulmod_shoup, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        const ulong d = n_randtest_not_zero(state) / 2 + 1;  // 0 < d < 2**(FLINT_BITS-1)
        const ulong a = n_randtest(state) % d;  // a must be < d
        const ulong b = n_randtest(state);  // b is arbitrary

        // mulmod_shoup
        const ulong a_pr = n_mulmod_precomp_shoup(a, d);
        const ulong r1 = n_mulmod_shoup(a, b, a_pr, d);

        // trivial mulmod
        ulong r2, q, p1, p2;
        umul_ppmm(p1, p2, a, b);
        p1 %= d;
        udiv_qrnnd(q, r2, p1, p2, d);

        result = (r1 == r2);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "a = %wu, b = %wu, d = %wu, a_pr = %wu\n"
                    "q = %wu, r1 = %wu, r2 = %wu\n",
                    a, b, d, a_pr, q, r1, r2);
    }

    TEST_FUNCTION_END(state);
}
