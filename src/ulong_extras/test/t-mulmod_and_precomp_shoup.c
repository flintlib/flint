/*
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_mulmod_and_precomp_shoup, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        const ulong d = n_randtest_not_zero(state) / 2 + 1;  // 0 < d < 2**(FLINT_BITS-1)
        const ulong a = n_randtest(state) % d;  // a must be < d
        const ulong b = n_randtest(state) % d;  // b must be < d

        // mulmod_and_precomp_shoup
        ulong a_pr_hi, a_pr_lo, ab, ab_precomp;

        n_mulmod_precomp_shoup_hi_lo(&a_pr_hi, &a_pr_lo, a, d);
        const ulong b_precomp = n_mulmod_precomp_shoup(b, d);
        n_mulmod_and_precomp_shoup(&ab, &ab_precomp, a, b, a_pr_hi, a_pr_lo, b_precomp, d);

        // check ab == a*b % n
        ulong r2, q, p1, p2;
        umul_ppmm(p1, p2, a, b);
        p1 %= d;
        udiv_qrnnd(q, r2, p1, p2, d);
        result = (ab == r2);

        // check precomp values:
        result = result & (a_pr_hi == n_mulmod_precomp_shoup(a, d));
        result = result & (a_pr_lo + a_pr_hi * d == 0L);
        result = result & (ab_precomp == n_mulmod_precomp_shoup(ab, d));

        if (!result)
            TEST_FUNCTION_FAIL("a = %wu, b = %wu, d = %wu, \n",
                    "a_pr_hi = %wu, a_pr_lo = %wu, b_precomp = %wu\n"
                    "q = %wu, r1 = %wu, r2 = %wu, ab_precomp = %wu\n",
                    a, b, d, a_pr_hi, a_pr_lo, b_precomp, q, ab, r2, ab_precomp);
    }

    TEST_FUNCTION_END(state);
}
