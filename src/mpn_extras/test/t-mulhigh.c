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

#define N_MAX (FLINT_MPN_MULHIGH_MUL_CUTOFF + 100)

#define rpcH (rpc + n - 1)

/* Defined in t-mulhigh.c and t-sqrhigh.c */
#ifndef lower_bound
# define lower_bound lower_bound
static ulong lower_bound(ulong n)
{
    /* These are calculated by hand lower bound for the returned limb */
    if (n < 3)
        return 0;
    else
        return 4 * n - 8;
}
#endif

TEST_FUNCTION_START(flint_mpn_mulhigh, state)
{
    slong ix;
    int result;

    mp_ptr rp, rpc, xp, yp;

    rp = flint_malloc(sizeof(mp_limb_t) * N_MAX);
    rpc = flint_malloc(2 * sizeof(mp_limb_t) * N_MAX);
    xp = flint_malloc(sizeof(mp_limb_t) * N_MAX);
    yp = flint_malloc(sizeof(mp_limb_t) * N_MAX);

    for (ix = 0; ix < 100000 * flint_test_multiplier(); ix++)
    {
        mp_limb_t borrow;
        mp_size_t n;

        /* Trigger full multiplication in mulhigh */
        if (n_randint(state, 1000) == 0)
            n = FLINT_MPN_MULHIGH_MUL_CUTOFF + n_randint(state, 50);
        else if (n_randint(state, 100) == 0)
            n = 1 + n_randint(state, FLINT_MPN_MULHIGH_MUL_CUTOFF);
        else
            n = 1 + n_randint(state, 2 * FLINT_MPN_MULHIGH_MULDERS_CUTOFF);

        mpn_random2(xp, n);
        mpn_random2(yp, n);

        rp[0] = flint_mpn_mulhigh_n(rp + 1, xp, yp, n);

        /* Check upper bound */
        flint_mpn_mul_n(rpc, xp, yp, n);

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
                    ix, n, xp, n, yp, n, rpcH, n + 1, rp, n + 1);

        /* Check lower bound */
        borrow = mpn_sub_1(rpcH, rpcH, n + 1, lower_bound(n));
        if (borrow)
            mpn_zero(rpcH, n + 1);

        result = (mpn_cmp(rp, rpcH, n + 1) >= 0);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "Smaller than lower bound!\n"
                    "ix = %wd\n"
                    "n = %wd\n"
                    "xp = %{ulong*}\n"
                    "yp = %{ulong*}\n"
                    "Lower bound: %{ulong*}\n"
                    "Got:         %{ulong*}\n",
                    ix, n, xp, n, yp, n, rpcH, n + 1, rp, n + 1);
    }

    flint_free(rp);
    flint_free(rpc);
    flint_free(xp);
    flint_free(yp);

    TEST_FUNCTION_END(state);
}
# undef N_MIN
# undef N_MAX
# undef N_MAX2
#else
TEST_FUNCTION_START(flint_mpn_mulhigh, state)
{
    TEST_FUNCTION_END_SKIPPED(state);
}
#endif
