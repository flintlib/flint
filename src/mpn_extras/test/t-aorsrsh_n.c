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

#define N_MIN                                        1
#define N_MAX   (FLINT_MPN_AORSRSH_FUNC_TAB_WIDTH -  1)
#define N_STOR  (FLINT_MPN_AORSRSH_FUNC_TAB_WIDTH + 10)

static
mp_limb_t mpn_addrsh_n(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n, unsigned int cnt)
{
    mpn_rshift(rp, yp, n, cnt);
    return mpn_add_n(rp, rp, xp, n);
}

static
mp_limb_t mpn_subrsh_n(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n, unsigned int cnt)
{
    mpn_rshift(rp, yp, n, cnt);
    return mpn_sub_n(rp, xp, rp, n);
}

TEST_FUNCTION_START(flint_mpn_aorsrsh_n, state)
{
#if FLINT_USE_AORSRSH_FUNC_TAB
    slong ix;

    for (ix = 0; ix < 10000 * flint_test_multiplier(); ix++)
    {
        int result;
        int type;
        int aliasing;
        unsigned int cnt;
        mp_limb_t cf, cg;
        mp_size_t n;
        mp_ptr fp, gp, xp, yp;

        n = N_MIN + n_randint(state, N_MAX - N_MIN + 1);
        if (n_randint(state, 1 << 10) == UWORD(0))
            n += N_STOR;

        /* 0: No aliasing
         * 1: fp = xp
         * 2: fp = yp */
        aliasing = n_randint(state, 3);

        fp = flint_malloc(sizeof(mp_limb_t) * n);
        gp = flint_malloc(sizeof(mp_limb_t) * n);
        xp = flint_malloc(sizeof(mp_limb_t) * n);
        yp = flint_malloc(sizeof(mp_limb_t) * n);

        flint_mpn_rrandom(xp, state, n);
        flint_mpn_rrandom(yp, state, n);
        cnt = 1 + n_randint(state, FLINT_BITS - 1);

        type = n_randint(state, 2);

        /* FIXME */
        if (n > N_MAX && aliasing == 1)
            aliasing = 0;

        if (type == 0)
        {
            if (aliasing == 0)
                cf = flint_mpn_addrsh_n(fp, xp, yp, n, cnt);
            else if (aliasing == 1)
            {
                flint_mpn_copyi(fp, xp, n);
                cf = flint_mpn_addrsh_n(fp, fp, yp, n, cnt);
            }
            else
            {
                flint_mpn_copyi(fp, yp, n);
                cf = flint_mpn_addrsh_n(fp, xp, fp, n, cnt);
            }
            cg = mpn_addrsh_n(gp, xp, yp, n, cnt);
        }
        else
        {
            if (aliasing == 0)
                cf = flint_mpn_subrsh_n(fp, xp, yp, n, cnt);
            else if (aliasing == 1)
            {
                flint_mpn_copyi(fp, xp, n);
                cf = flint_mpn_subrsh_n(fp, fp, yp, n, cnt);
            }
            else
            {
                flint_mpn_copyi(fp, yp, n);
                cf = flint_mpn_subrsh_n(fp, xp, fp, n, cnt);
            }
            cg = mpn_subrsh_n(gp, xp, yp, n, cnt);
        }

        result = (cf == cg && mpn_cmp(fp, gp, n) == 0);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "function: %s\n"
                    "aliasing: %s\n"
                    "ix = %wd\n"
                    "n = %wd\n"
                    "cnt = %u\n"
                    "xp = %{ulong*}\n"
                    "yp = %{ulong*}\n"
                    "FLINT (cy = %wu): %{ulong*}\n"
                    "GMP   (cy = %wu): %{ulong*}\n",
                    type == 0 ? "flint_mpn_addrsh_n" : "flint_mpn_subrsh_n",
                    aliasing == 0 ? "none" : (aliasing == 1 ? "rp = xp" : "rp = yp"),
                    ix, n, cnt, xp, n, yp, n, cf, fp, n, cg, gp, n);

        flint_free(fp);
        flint_free(gp);
        flint_free(xp);
        flint_free(yp);
    }

    TEST_FUNCTION_END(state);
#else
    TEST_FUNCTION_END_SKIPPED(state);
#endif
}
#undef N_MIN
#undef N_MAX
#undef N_STOR
