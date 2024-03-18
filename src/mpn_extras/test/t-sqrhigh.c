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

/* TODO: Remove this preprocessor conditional */
#if FLINT_HAVE_NATIVE_mpn_mulhigh_basecase

#define N_MAX FLINT_MAX(FLINT_MPN_SQRHIGH_SQR_CUTOFF + 50, 2 * FLINT_MPN_SQRHIGH_MULDERS_CUTOFF)

#define rpcH (rpc + n - 1)

TEST_FUNCTION_START(flint_mpn_sqrhigh, state)
{
    slong ix;
    int result;

    mp_ptr rp, rpc, xp;

    rp = flint_malloc(sizeof(mp_limb_t) * (N_MAX + 1));
    rpc = flint_malloc(2 * sizeof(mp_limb_t) * N_MAX);
    xp = flint_malloc(sizeof(mp_limb_t) * N_MAX);

    for (ix = 0; ix < 100000 * flint_test_multiplier(); ix++)
    {
        mp_limb_t borrow;
        mp_size_t n;
        mp_limb_t lb;

        /* Trigger full multiplication in mulhigh */
        if (n_randint(state, 1000) == 0)
            n = 1 + FLINT_MPN_SQRHIGH_SQR_CUTOFF + n_randint(state, 50);
        else if (n_randint(state, 100) == 0)
            n = 1 + n_randint(state, FLINT_MPN_SQRHIGH_SQR_CUTOFF);
        else
            n = 1 + n_randint(state, 2 * FLINT_MPN_SQRHIGH_MULDERS_CUTOFF);

        mpn_random2(xp, n);

        rp[0] = flint_mpn_sqrhigh(rp + 1, xp, n);

        /* Check upper bound */
        flint_mpn_sqr(rpc, xp, n);

        result = (mpn_cmp(rp, rpcH, n + 1) <= 0);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "Bigger than upper bound!\n"
                    "ix = %wd\n"
                    "n = %wd\n"
                    "xp = %{ulong*}\n"
                    "yp = %{ulong*}\n"
                    "Upper bound: %{ulong*}\n"
                    "Got:         %{ulong*}\n",
                    ix, n, xp, n, rpcH, n + 1, rp, n + 1);

        /* Check lower bound */
        lb = 2 * n;

        borrow = mpn_sub_1(rpcH, rpcH, n + 1, lb);
        if (borrow)
            mpn_zero(rpcH, n + 1);

        result = (mpn_cmp(rp, rpcH, n + 1) >= 0);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "Smaller than lower bound!\n"
                    "ix = %wd\n"
                    "n = %wd\n"
                    "xp = %{ulong*}\n"
                    "Lower bound: %{ulong*}\n"
                    "Got:         %{ulong*}\n",
                    ix, n, xp, n, rpcH, n + 1, rp, n + 1);
    }

    flint_free(rp);
    flint_free(rpc);
    flint_free(xp);

    TEST_FUNCTION_END(state);
}
# undef N_MAX
# undef rpcH
#else
TEST_FUNCTION_START(flint_mpn_sqrhigh, state)
{
    TEST_FUNCTION_END_SKIPPED(state);
}
#endif
