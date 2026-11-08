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

#define N_MIN 1

#define WANT_STOR 1

#if WANT_STOR
# define MAX_ALLOC_SIZE (FLINT_MPN_SQRHIGH_SQR_CUTOFF + 50)
#else
# define N_MAX FLINT_MPN_SQRHIGH_FUNC_TAB_WIDTH
# define MAX_ALLOC_SIZE (2 * N_MAX)
#endif

#define GET_N_FULL_SQR(state) (N_MIN + FLINT_MPN_SQRHIGH_SQR_CUTOFF + n_randint(state, 50 - N_MIN + 1))
#define GET_N_LARGE(state) (N_MIN + n_randint(state, FLINT_MPN_SQRHIGH_SQR_CUTOFF - N_MIN + 1))
#define GET_N_SMALL(state) (N_MIN + n_randint(state, 2 * FLINT_MPN_SQRHIGH_MULDERS_CUTOFF - N_MIN + 1))

#define rpcH (rpc + n - 1)


/* Defined in t-mulhigh_n.c and t-sqrhigh.c */
#ifndef lower_bound
# define lower_bound(n) FLINT_MPN_MULHIGH_N_GUARD_ERR_BOUND(n)

/* Inputs maximising the dropped low half-products on the diagonal
   p + q = n - 2 while keeping the dropped diagonal p + q = n - 3 nearly
   maximal: the deficit reaches 2n - 4 units of the guard limb, the
   maximum possible. For
   squaring the two values alternate, which pairs them on the diagonal
   when n is odd. */
static void adversarial(mp_ptr xp, mp_ptr yp, mp_size_t n, flint_rand_t state)
{
    mp_size_t i;
    mp_limb_t a = -(UWORD(1) << (FLINT_BITS / 2)) - 1, b = -(UWORD(1) << (FLINT_BITS / 2)) + 1;

    for (i = 0; i < n; i++)
    {
        if (yp == xp)
            xp[i] = (i & 1) ? a : b;
        else
            xp[i] = a, yp[i] = b;
    }

    /* the top limbs do not matter for the dropped part; randomising them
       exercises the carry propagation */
    if (n_randint(state, 2))
    {
        xp[n - 1] = n_randtest(state) | (UWORD(1) << (FLINT_BITS - 1));
        if (yp != xp)
            yp[n - 1] = n_randtest(state) | (UWORD(1) << (FLINT_BITS - 1));
    }
}
#endif

TEST_FUNCTION_START(flint_mpn_sqrhigh_n, state)
{
    slong ix;
    int result;
    mp_ptr rp, rpc, xp;

    rpc = flint_malloc(2 * sizeof(mp_limb_t) * MAX_ALLOC_SIZE);

    for (ix = 0; ix < 10000 * flint_test_multiplier(); ix++)
    {
        mp_limb_t borrow;
        mp_size_t n;

#if WANT_STOR
        /* Trigger full multiplication in mulhigh */
        if (n_randint(state, 1000) == 0)
            n = GET_N_FULL_SQR(state);
        else if (n_randint(state, 100) == 0)
            n = GET_N_LARGE(state);
        else
            n = GET_N_SMALL(state);
#else
        n = N_MIN + n_randint(state, N_MAX - N_MIN + 1);
#endif

        xp = flint_malloc(sizeof(mp_limb_t) * n);

        if (n_randint(state, 4) == 0)
            adversarial(xp, xp, n, state);
        else
            flint_mpn_rrandom(xp, state, n);

        if (n <= FLINT_MAX(FLINT_MPN_SQRHIGH_MULDERS_CUTOFF, FLINT_MPN_SQRHIGH_FUNC_TAB_WIDTH))
        {
            /* Check that it is coherent with generic version */
            mp_limb_t ret0, ret1;

            rp = flint_malloc(sizeof(mp_limb_t) * n);

            ret0 = _flint_mpn_mulhigh_n_naive(rpc, xp, xp, n);
            ret1 = flint_mpn_sqrhigh_n(rp, xp, n);

            result = (mpn_cmp(rp, rpc, n) == 0 && ret0 == ret1);
            if (!result)
                TEST_FUNCTION_FAIL(
                        "Basecase not coherent with generic version.\n"
                        "ix = %wd\n"
                        "n = %wd\n"
                        "xp = %{ulong*}\n"
                        "Expected ret: %{ulong}\n"
                        "Got ret:      %{ulong}\n"
                        "Expected limbs: %{ulong*}\n"
                        "Got limbs:      %{ulong*}\n",
                        ix, n, xp, n, ret0, ret1, rpc, n, rp, n);
            flint_free(rp);
        }

        {
            /* Check that it is coherent with bounds */
            rp = flint_malloc(sizeof(mp_limb_t) * (n + 1));

            rp[0] = flint_mpn_sqrhigh_n(rp + 1, xp, n);
            flint_mpn_sqr(rpc, xp, n);

            result = (mpn_cmp(rp, rpcH, n + 1) <= 0);
            if (!result)
                TEST_FUNCTION_FAIL(
                        "Bigger than upper bound!\n"
                        "ix = %wd\n"
                        "n = %wd\n"
                        "xp = %{ulong*}\n"
                        "Upper bound: %{ulong*}\n"
                        "Got:         %{ulong*}\n",
                        ix, n, xp, n, rpcH, n + 1, rp, n + 1);

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
                        "Lower bound: %{ulong*}\n"
                        "Got:         %{ulong*}\n",
                        ix, n, xp, n, rpcH, n + 1, rp, n + 1);
        }

        /* aliasing: the output may be the operand */
        {
            mp_ptr xc = flint_malloc(sizeof(mp_limb_t) * n);
            mp_limb_t ret1;

            flint_mpn_copyi(xc, xp, n);
            ret1 = flint_mpn_sqrhigh_n(xc, xc, n);
            result = (mpn_cmp(xc, rp + 1, n) == 0 && ret1 == rp[0]);
            if (!result)
                TEST_FUNCTION_FAIL("aliasing res = x failed\nn = %wd\n", n);

            flint_free(xc);
        }

        flint_free(xp);
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

#undef GET_N_FULL_SQR
#undef GET_N_LARGE
#undef GET_N_SMALL
