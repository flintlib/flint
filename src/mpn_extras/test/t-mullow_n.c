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
#define N_MAX 30

TEST_FUNCTION_START(flint_mpn_mullow_n, state)
{
    slong ix;
    int result;

    mp_ptr rp, rpf, xp, yp;

    for (ix = 0; ix < 100000 * flint_test_multiplier(); ix++)
    {
        mp_limb_t ret;
        mp_size_t n;

        n = N_MIN + n_randint(state, N_MAX - N_MIN + 1);

        rp = flint_malloc(sizeof(mp_limb_t) * n);
        rpf = flint_malloc(2 * sizeof(mp_limb_t) * n);
        xp = flint_malloc(sizeof(mp_limb_t) * n);
        yp = flint_malloc(sizeof(mp_limb_t) * n);

        flint_mpn_rrandom(xp, state, n);
        flint_mpn_rrandom(yp, state, n);

        ret = flint_mpn_mullow_n(rp, xp, yp, n);
        flint_mpn_mul_n(rpf, xp, yp, n);

        result = (mpn_cmp(rp, rpf, n) == 0 && ret == rpf[n]);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "ix = %wd\n"
                    "n = %wd\n"
                    "xp = %{ulong*}\n"
                    "yp = %{ulong*}\n"
                    "Exp ret: %{ulong}\n"
                    "Got ret: %{ulong}\n"
                    "Expected: %{ulong*}\n"
                    "Got:      %{ulong*}\n",
                    ix, n, xp, n, yp, n, rpf[n], ret, rpf, n, rp, n);

        flint_free(rp);
        flint_free(rpf);
        flint_free(xp);
        flint_free(yp);
    }

    TEST_FUNCTION_END(state);
}

#undef N_MIN
#undef N_MAX
