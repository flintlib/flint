/*
    Copyright (C) 2024 Albin Ahlbäck
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"

#define N_MIN 1

#define WANT_STOR 1

#if WANT_STOR
# define MAX_ALLOC_SIZE (FLINT_MPN_MULHIGH_MUL_CUTOFF + 50)
#else
# define N_MAX FLINT_MPN_MULHIGH_FUNC_TAB_WIDTH
# define MAX_ALLOC_SIZE (2 * N_MAX)
#endif

#define GET_N_FULL_MUL(state) (N_MIN + FLINT_MPN_MULHIGH_MUL_CUTOFF + n_randint(state, 50 - N_MIN + 1))
#define GET_N_LARGE(state) (N_MIN + n_randint(state, FLINT_MPN_MULHIGH_MUL_CUTOFF - N_MIN + 1))
#define GET_N_SMALL(state) (N_MIN + n_randint(state, 2 * FLINT_MPN_MULHIGH_MULDERS_CUTOFF - N_MIN + 1))

#define rpcH (rpc + n - 1)

/* Defined in t-mulhigh.c and t-sqrhigh.c */
#ifndef generic_mulhigh_basecase
# define generic_mulhigh_basecase generic_mulhigh_basecase
static mp_limb_t
generic_mulhigh_basecase(mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n)
{
    mp_limb_t low;
    mp_limb_t a, b;

    if (n == 1)
    {
        umul_ppmm(rp[0], low, up[0], vp[0]);
    }
    else if (n == 2)
    {
        FLINT_MPN_MUL_2X2(rp[1], rp[0], low, b, up[1], up[0], vp[1], vp[0]);
    }
    else if (n == 3)
    {
        NN_DOTREV_S3_1X1_HIGH(b, a, up, vp, 2);
        NN_DOTREV_S3_A3_1X1(b, a, low, 0, b, a, up, vp, 3);
        NN_DOTREV_S3_A3_1X1(b, a, rp[0], 0, b, a, up + 1, vp + 1, 2);
        NN_ADDMUL_S2_A2_1X1(rp[2], rp[1], b, a, up[2], vp[2]);
    }
    else
    {
        mp_limb_t tp[2];
        slong ix;

        NN_DOTREV_S3_1X1_HIGH(b, a, up, vp, n - 1);
        NN_DOTREV_S3_A3_1X1(b, a, low, 0, b, a, up, vp, n);
        tp[0] = a;
        tp[1] = b;

        umul_ppmm(rp[1], rp[0], up[n - 1], vp[1]);
        for (ix = 2; ix < n; ix++)
            rp[ix] = mpn_addmul_1(rp, up + n - ix, ix, vp[ix]);

        mpn_add(rp, rp, n, tp, 2);
    }

    return low;
}
#endif

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

TEST_FUNCTION_START(flint_mpn_mulhigh_n, state)
{
    slong ix;
    int result;
    mp_ptr rp, rpc, xp, yp;

    rpc = flint_malloc(2 * sizeof(mp_limb_t) * MAX_ALLOC_SIZE);

    for (ix = 0; ix < 100000 * flint_test_multiplier(); ix++)
    {
        mp_limb_t borrow;
        mp_size_t n;

#if WANT_STOR
        /* Trigger full multiplication in mulhigh */
        if (n_randint(state, 1000) == 0)
            n = GET_N_FULL_MUL(state);
        else if (n_randint(state, 100) == 0)
            n = GET_N_LARGE(state);
        else
            n = GET_N_SMALL(state);
#else
        n = N_MIN + n_randint(state, N_MAX - N_MIN + 1);
#endif

        xp = flint_malloc(sizeof(mp_limb_t) * n);
        yp = flint_malloc(sizeof(mp_limb_t) * n);

        flint_mpn_rrandom(xp, state, n);
        flint_mpn_rrandom(yp, state, n);

        if (n <= FLINT_MAX(FLINT_MPN_MULHIGH_MULDERS_CUTOFF, FLINT_MPN_MULHIGH_FUNC_TAB_WIDTH))
        {
            /* Check that it is coherent with generic version */
            mp_limb_t ret0, ret1;

            rp = flint_malloc(sizeof(mp_limb_t) * n);

            ret0 = generic_mulhigh_basecase(rpc, xp, yp, n);
            ret1 = flint_mpn_mulhigh_n(rp, xp, yp, n);

            result = (mpn_cmp(rp, rpc, n) == 0 && ret0 == ret1);
            if (!result)
                TEST_FUNCTION_FAIL(
                        "Basecase not coherent with generic version.\n"
                        "ix = %wd\n"
                        "n = %wd\n"
                        "xp = %{ulong*}\n"
                        "yp = %{ulong*}\n"
                        "Expected ret: %{ulong}\n"
                        "Got ret:      %{ulong}\n"
                        "Expected limbs: %{ulong*}\n"
                        "Got limbs:      %{ulong*}\n",
                        ix, n, xp, n, yp, n, ret0, ret1, rpc, n, rp, n);
        }
        else
        {
            /* Check that it is coherent with bounds */
            rp = flint_malloc(sizeof(mp_limb_t) * (n + 1));

            rp[0] = flint_mpn_mulhigh_n(rp + 1, xp, yp, n);
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

        flint_free(xp);
        flint_free(yp);
        flint_free(rp);
    }

    flint_free(rpc);

    TEST_FUNCTION_END(state);
}

#undef N_MIN
#undef N_MAX
#undef WANT_STOR
#undef rpcH
#undef MAX_ALLOC_SIZE

#undef GET_N_FULL_MUL
#undef GET_N_LARGE
#undef GET_N_SMALL
