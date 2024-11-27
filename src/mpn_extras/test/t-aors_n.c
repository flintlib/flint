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

#define N_MIN                                     1
#define N_MAX   (FLINT_MPN_AORS_FUNC_TAB_WIDTH -  1)
#define N_STOR  (FLINT_MPN_AORS_FUNC_TAB_WIDTH + 10)

TEST_FUNCTION_START(flint_mpn_aors_n, state)
{
#if FLINT_USE_AORS_FUNC_TAB
    slong ix;

    for (ix = 0; ix < 10000 * flint_test_multiplier(); ix++)
    {
        int result;
        int type;
        mp_limb_t cf, cg;
        mp_size_t n;
        mp_ptr fp, gp, xp, yp;

        n = N_MIN + n_randint(state, N_MAX - N_MIN + 1);
        if (n_randint(state, 1 << 10) == UWORD(0))
            n += N_STOR;

        fp = flint_malloc(sizeof(mp_limb_t) * n);
        gp = flint_malloc(sizeof(mp_limb_t) * n);
        xp = flint_malloc(sizeof(mp_limb_t) * n);
        yp = flint_malloc(sizeof(mp_limb_t) * n);

        flint_mpn_rrandom(xp, state, n);
        flint_mpn_rrandom(yp, state, n);

        type = n_randint(state, 2);

        if (type == 0)
        {
            cf = flint_mpn_add_n(fp, xp, yp, n);
            cg = mpn_add_n(gp, xp, yp, n);
        }
        else
        {
            cf = flint_mpn_sub_n(fp, xp, yp, n);
            cg = mpn_sub_n(gp, xp, yp, n);
        }

        result = (cf == cg && mpn_cmp(fp, gp, n) == 0);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "%s:\n"
                    "ix = %wd\n"
                    "n = %wd\n"
                    "xp = %{ulong*}\n"
                    "yp = %{ulong*}\n"
                    "FLINT (cy = %wu): %{ulong*}\n"
                    "GMP   (cy = %wu): %{ulong*}\n",
                    type == 0 ? "flint_mpn_add_n" : "flint_mpn_sub_n",
                    ix, n, xp, n, yp, n, cf, fp, n, cg, gp, n + 1);

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
