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

    /* division-free precomputations (preinv) against the ones by division,
       for any modulus 1 <= n < 2^FLINT_BITS and a < n */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        ulong n, ninv, a, q1, q2, q3, q4, r1, r2;

        switch (n_randint(state, 4))
        {
            case 0:
                n = n_randtest_not_zero(state);
                break;
            case 1:   /* 2^k, 2^k - 1, 2^k + 1 */
                n = UWORD(1) << n_randint(state, FLINT_BITS);
                n = n + (n_randint(state, 3) - 1);
                if (n == 0)
                    n = 1;
                break;
            case 2:
                n = n_randlimb(state) | (UWORD(1) << (FLINT_BITS - 1));
                break;
            default:
                n = n_randbits(state, 1 + n_randint(state, FLINT_BITS));
                break;
        }

        switch (n_randint(state, 4))
        {
            case 0: a = 0; break;
            case 1: a = n - 1; break;
            case 2: a = n_randint(state, FLINT_MIN(n, 3)); break;
            default: a = n_randint(state, n); break;
        }

        ninv = n_preinvert_limb(n);

        q1 = n_mulmod_precomp_shoup(a, n);
        q2 = n_mulmod_precomp_shoup_preinv(a, n, ninv);
        n_mulmod_precomp_shoup_quo_rem(&q3, &r1, a, n);
        n_mulmod_precomp_shoup_quo_rem_preinv(&q4, &r2, a, n, ninv);

        result = (q2 == q1 && q3 == q1 && q4 == q1 && r2 == r1);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "n = %wu, a = %wu\n"
                    "division: quo = %wu, %wu, rem = %wu\n"
                    "preinv: quo = %wu, %wu, rem = %wu\n",
                    n, a, q1, q3, r1, q2, q4, r2);
    }

    TEST_FUNCTION_END(state);
}
