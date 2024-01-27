/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"

#if FLINT_HAVE_MPN_MULHIGH_BASECASE

#define N_MIN 6
#define N_MAX 64

TEST_FUNCTION_START(flint_mpn_mulhigh_basecase, state)
{
    slong ix;
    int result;

    for (ix = 0; ix < 100000 * flint_test_multiplier(); ix++)
    {
        mp_limb_t rp[N_MAX + 1];
        mp_limb_t rp_upperbound[2 * N_MAX];
        mp_limb_t rp_lowerbound[2 * N_MAX];
        mp_limb_t borrow;
        mp_limb_t xp[N_MAX];
        mp_limb_t yp[N_MAX];
        mp_size_t n;
        mp_limb_t lb;

        n = N_MIN + n_randint(state, N_MAX - N_MIN + 1);

        mpn_random2(xp, N_MAX);
        mpn_random2(yp, N_MAX);

        rp[0] = flint_mpn_mulhigh_basecase(rp + 1, xp, yp, n);

        /* Check upper bound */
        flint_mpn_mul_n(rp_upperbound, xp, yp, n);

        result = (mpn_cmp(rp, rp_upperbound + n - 1, n + 1) <= 0);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "rp > rp_upperbound\n"
                    "ix = %wd\n"
                    "n = %wd\n"
                    "xp = %{ulong*}\n"
                    "yp = %{ulong*}\n"
                    "rp            = %{ulong*}\n"
                    "rp_upperbound = %{ulong*}\n",
                    ix, n, xp, n, yp, n, rp, n + 1, rp_upperbound + n - 1, n + 1);

        /* Check lower bound */
        lb = n + 2;
        memcpy(rp_lowerbound, rp_upperbound, 2 * n * sizeof(mp_limb_t));

        borrow = mpn_sub_1(rp_lowerbound + n - 1, rp_lowerbound + n - 1, n + 1, lb);
        if (borrow)
            mpn_zero(rp_lowerbound, 2 * n);

        result = (mpn_cmp(rp, rp_lowerbound + n - 1, n + 1) >= 0);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "rp < rp_lowerbound\n"
                    "ix = %wd\n"
                    "n = %wd\n"
                    "xp = %{ulong*}\n"
                    "yp = %{ulong*}\n"
                    "rp            = %{ulong*}\n"
                    "rp_lowerbound = %{ulong*}\n",
                    ix, n, xp, n, yp, n, rp, n + 1, rp_lowerbound + n - 1, n + 1);
    }

#undef N_MAX
#define N_MAX 12

    /* Check against hardcoded versions */
    for (ix = 0; ix < 10000 * flint_test_multiplier(); ix++)
    {
        mp_limb_t rp1[N_MAX + 1];
        mp_limb_t rp2[N_MAX + 1];
        mp_limb_t xp[N_MAX + 10];
        mp_limb_t yp[N_MAX + 10];
        mp_size_t n;

        n = N_MIN + n_randint(state, N_MAX - N_MIN + 1);

        mpn_random2(xp, N_MAX + 10);
        mpn_random2(yp, N_MAX + 10);

        rp1[0] = flint_mpn_mulhigh_n(rp1 + 1, xp, yp, n);
        rp2[0] = flint_mpn_mulhigh_basecase(rp2 + 1, xp, yp, n);

        result = (rp1[0] == rp2[0]);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "Return values are not the same\n"
                    "ix = %wd\n"
                    "n = %wd\n"
                    "xp = %{ulong*}\n"
                    "yp = %{ulong*}\n"
                    "Expected: %{ulong*}\n"
                    "Got:      %{ulong*}\n",
                    ix, n, xp, n, yp, n, rp1, n + 1, rp2, n + 1);

        result = (mpn_cmp(rp1, rp2, n + 1) == 0);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "Result not the same\n"
                    "ix = %wd\n"
                    "n = %wd\n"
                    "xp = %{ulong*}\n"
                    "yp = %{ulong*}\n"
                    "Expected: %{ulong*}\n"
                    "Got:      %{ulong*}\n",
                    ix, n, xp, n, yp, n, rp1, n + 1, rp2, n + 1);
    }

    TEST_FUNCTION_END(state);
}
# undef N_MIN
# undef N_MAX
#else
TEST_FUNCTION_START(flint_mpn_mulhigh_basecase, state)
{
    TEST_FUNCTION_END_SKIPPED(state);
}
#endif
