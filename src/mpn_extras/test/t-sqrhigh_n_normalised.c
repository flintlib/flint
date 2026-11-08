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
#define N_MAX 64

TEST_FUNCTION_START(flint_mpn_sqrhigh_n_normalised, state)
{
    slong ix;
    int result;

    for (ix = 0; ix < 10000 * flint_test_multiplier(); ix++)
    {
        mp_ptr rp_n, rp_u, xp;
        mp_size_t n;
        mp_limb_pair_t res_norm;
        mp_limb_t retlimb, normalised;

        n = N_MIN + n_randint(state, N_MAX - N_MIN + 1);

        rp_n = flint_malloc(sizeof(mp_limb_t) * (n + 1));
        rp_u = flint_malloc(sizeof(mp_limb_t) * (n + 1));
        xp = flint_malloc(sizeof(mp_limb_t) * n);

        if (n_randint(state, 4) == 0)
        {
            /* adversarial: deficit close to the bound (see t-sqrhigh.c) */
            mp_size_t i;
            for (i = 0; i < n; i++)
                xp[i] = (i & 1) ? -(UWORD(1) << (FLINT_BITS / 2)) - 1 : -(UWORD(1) << (FLINT_BITS / 2)) + 1;
        }
        else
            flint_mpn_rrandom(xp, state, n);
        xp[n - 1] |= (UWORD(1) << (FLINT_BITS - 1));

        rp_u[0] = flint_mpn_sqrhigh_n(rp_u + 1, xp, n);
        res_norm = flint_mpn_sqrhigh_n_normalised(rp_n + 1, xp, n);

        retlimb = res_norm.m1;
        normalised = res_norm.m2;
        rp_n[0] = retlimb;

        result = ((rp_n[n] & (UWORD(1) << (FLINT_BITS - 1))) != UWORD(0));
        if (!result)
            TEST_FUNCTION_FAIL(
                    "Top bit not set in normalised result\n"
                    "ix = %wd\n"
                    "n = %wd\n"
                    "xp = %{ulong*}\n"
                    "rp_n = %{ulong*}\n"
                    "rp_u = %{ulong*}\n",
                    ix, n, xp, n, rp_n, n + 1, rp_u, n + 1);

        if (normalised)
        {
            result = (mpn_lshift(rp_u, rp_u, n + 1, 1) == 0);
            result = result && (mpn_cmp(rp_n, rp_u, n + 1) == 0);
            if (!result)
                TEST_FUNCTION_FAIL(
                        "rp_n != rp_u << 1 when normalised\n"
                        "ix = %wd\n"
                        "n = %wd\n"
                        "xp = %{ulong*}\n"
                        "rp_n = %{ulong*}\n"
                        "rp_u = %{ulong*}\n",
                        ix, n, xp, n, rp_n, n + 1, rp_u, n + 1);
        }
        else
        {
            result = (mpn_cmp(rp_n, rp_u, n + 1) == 0);
            if (!result)
                TEST_FUNCTION_FAIL(
                        "rp_n != rp_u when unnormalised\n"
                        "ix = %wd\n"
                        "n = %wd\n"
                        "xp = %{ulong*}\n"
                        "rp_n = %{ulong*}\n"
                        "rp_u = %{ulong*}\n",
                        ix, n, xp, n, rp_n, n + 1, rp_u, n + 1);
        }


        /* Check the deficit against the exact (normalised) product */
        {
            mp_ptr full = flint_malloc(sizeof(mp_limb_t) * 2 * n);
            mp_limb_t bound = FLINT_MPN_SQRHIGH_N_NORMALISED_GUARD_ERR_BOUND(n);
            mp_size_t i;

            flint_mpn_sqr(full, xp, n);
            if (normalised)
                mpn_lshift(full, full, 2 * n, 1);

            /* exact top n + 1 limbs minus computed: must be in [0, bound] */
            result = (mpn_sub_n(full + n - 1, full + n - 1, rp_n, n + 1) == 0);
            for (i = n; i < 2 * n; i++)
                result = result && (full[i] == 0);
            result = result && (full[n - 1] <= bound);
            if (!result)
                TEST_FUNCTION_FAIL(
                        "Deficit out of bounds\n"
                        "ix = %wd\n"
                        "n = %wd\n"
                        "xp = %{ulong*}\n"
                        "deficit = %{ulong}, bound = %{ulong}\n",
                        ix, n, xp, n, full[n - 1], bound);
            flint_free(full);
        }

        flint_free(rp_n);
        flint_free(rp_u);
        flint_free(xp);
    }

    TEST_FUNCTION_END(state);
}
#undef N_MIN
#undef N_MAX
