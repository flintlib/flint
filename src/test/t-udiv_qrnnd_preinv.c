/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "test_helpers.h"

TEST_FUNCTION_START(udiv_qrnnd_preinv, state)
{
    int i, result;

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        ulong d, dinv, nh, nl, q1, r1, q2, r2, norm;

        do
        {
            d = n_randtest_not_zero(state);
            nh = n_randtest(state);
            norm = flint_clz(d);
            d <<= norm;
        } while (nh >= d);
        nl = n_randtest(state);

        dinv = n_preinvert_limb_prenorm(d);

        udiv_qrnnd_preinv(q1, r1, nh, nl, d, dinv);
        udiv_qrnnd(q2, r2, nh, nl, d);

        result = ((q1 == q2) && (r1 == r2));
        if (!result)
            TEST_FUNCTION_FAIL(
                    "nh = %wu, nl = %wu, d = %wu, dinv = %wu\n"
                    "q1 = %wu, q2 = %wu, r1 = %wu, r2 = %wu\n",
                    nh, nl, d, dinv, q1, q2, r1, r2);
    }

    /*
        Edge cases: d = 2^(FLINT_BITS - 1), 2^(FLINT_BITS - 1) + 1, B - 2,
        B - 1 or random, nh in {0, 1, d - 2, d - 1} or random, nl in
        {0, 1, B - 1} or random; also with q and r overwriting nl and nh,
        as in a loop over the limbs of a dividend.
    */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        ulong d, dinv, nh, nl, q1, r1, q2, r2;

        switch (n_randint(state, 5))
        {
            case 0: d = UWORD(1) << (FLINT_BITS - 1); break;
            case 1: d = (UWORD(1) << (FLINT_BITS - 1)) + 1; break;
            case 2: d = UWORD_MAX - 1; break;
            case 3: d = UWORD_MAX; break;
            default: d = n_randlimb(state) | (UWORD(1) << (FLINT_BITS - 1));
        }

        switch (n_randint(state, 5))
        {
            case 0: nh = 0; break;
            case 1: nh = 1; break;
            case 2: nh = d - 2; break;
            case 3: nh = d - 1; break;
            default: nh = n_randlimb(state) % d;
        }

        switch (n_randint(state, 4))
        {
            case 0: nl = 0; break;
            case 1: nl = 1; break;
            case 2: nl = UWORD_MAX; break;
            default: nl = n_randlimb(state);
        }

        dinv = n_preinvert_limb_prenorm(d);

        udiv_qrnnd(q2, r2, nh, nl, d);
        udiv_qrnnd_preinv(q1, r1, nh, nl, d, dinv);

        if (q1 != q2 || r1 != r2)
            TEST_FUNCTION_FAIL(
                    "edge cases\n"
                    "nh = %wu, nl = %wu, d = %wu, dinv = %wu\n"
                    "q1 = %wu, q2 = %wu, r1 = %wu, r2 = %wu\n",
                    nh, nl, d, dinv, q1, q2, r1, r2);

        r1 = nh;
        q1 = nl;
        udiv_qrnnd_preinv(q1, r1, r1, q1, d, dinv);

        if (q1 != q2 || r1 != r2)
            TEST_FUNCTION_FAIL(
                    "aliasing\n"
                    "nh = %wu, nl = %wu, d = %wu, dinv = %wu\n"
                    "q1 = %wu, q2 = %wu, r1 = %wu, r2 = %wu\n",
                    nh, nl, d, dinv, q1, q2, r1, r2);
    }

    TEST_FUNCTION_END(state);
}
