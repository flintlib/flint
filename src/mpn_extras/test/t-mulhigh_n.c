/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "test_helpers.h"
#include "mpn_extras.h"

#if FLINT_MPN_MULHIGH_N_FUNC_TAB_WIDTH != 0

# define TEST_MPFR 0
# define N_MAX FLINT_MPN_MULHIGH_N_FUNC_TAB_WIDTH

void mpfr_mulhigh_n(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n);

TEST_FUNCTION_START(flint_mpn_mulhigh_n, state)
{
    slong ix;
    int result;

    for (ix = 0; ix < 100000 * flint_test_multiplier(); ix++)
    {
        mp_limb_t rp[N_MAX + 1] = {UWORD(0)};
        mp_limb_t rp_upperbound[2 * N_MAX] = {UWORD(0)};
        mp_limb_t rp_lowerbound[2 * N_MAX] = {UWORD(0)};
#if TEST_MPFR
        mp_limb_t rp_mpfr[2 * N_MAX] = {UWORD(0)};
#endif
        mp_limb_t borrow;
        mp_limb_t xp[N_MAX];
        mp_limb_t yp[N_MAX];
        mp_size_t n;
        mp_limb_t lb;

        n = 1 + n_randint(state, N_MAX);

        mpn_random2(xp, n);
        mpn_random2(yp, n);

        rp[0] = flint_mpn_mulhigh_n(rp + 1, xp, yp, n);

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
        memcpy(rp_lowerbound, rp_upperbound, 2 * n * sizeof(mp_limb_t));

        if (n <= 2)
            lb = 0; /* Multiplication should be "exact" */
        else if (n == 3)
            lb = 2;
        else if (n == 4)
            lb = 4;
        else if (n < 6)
            lb = n + 1;
        else
            lb = n + 2;

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

#if TEST_MPFR
        /* Check against MPFR */
        mpfr_mulhigh_n(rp_mpfr, xp, yp, n);

        result = (mpn_cmp(rp + n - 1, rp_mpfr + n - 1, n + 1) == 0);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "flint_mpn_mulhigh_n does not agree with mpfr_mulhigh_n\n"
                    "ix = %wd\n"
                    "n = %wd\n"
                    "xp = %{ulong*}\n"
                    "yp = %{ulong*}\n"
                    "rp      = %{ulong*}\n"
                    "rp_mpfr = %{ulong*}\n",
                    ix, n, xp, n, yp, n, rp + n - 1, n + 1, rp_mpfr + n - 1, n + 1);
#endif
    }

    TEST_FUNCTION_END(state);
}
#else
TEST_FUNCTION_START(flint_mpn_mulhigh_n, state)
{
    TEST_FUNCTION_END_SKIPPED(state);
}
#endif
