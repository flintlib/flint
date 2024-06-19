/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"

#define N_MIN 1
#define N_MAX 20

#define WANT_BIG 1
#define N_MIN_STOR (N_MAX + 1)
#define N_MAX_STOR 1000

TEST_FUNCTION_START(flint_mpn_mul, state)
{
    slong ix;

    for (ix = 0; ix < 100000 * flint_test_multiplier(); ix++)
    {
        slong m, n;
        mp_ptr xp, yp, rp1, rp2;
        mp_limb_t ret1, ret2;

        m = N_MIN + n_randint(state, N_MAX - N_MIN + 1);
#if WANT_BIG
        if (n_randint(state, 10000) == 0)
            m = N_MIN_STOR + n_randint(state, N_MAX_STOR - N_MIN_STOR + 1);
#endif
        n = 1 + n_randint(state, m);

        xp = flint_malloc(sizeof(mp_limb_t) * m);
        yp = flint_malloc(sizeof(mp_limb_t) * n);
        rp1 = flint_malloc(sizeof(mp_limb_t) * (m + n));
        rp2 = flint_malloc(sizeof(mp_limb_t) * (m + n));

        flint_mpn_rrandom(xp, state, m);
        flint_mpn_rrandom(yp, state, n);

        flint_mpn_rrandom(rp2, state, m + n);

        ret1 = mpn_mul(rp1, xp, m, yp, n);
        ret2 = flint_mpn_mul(rp2, xp, m, yp, n);

        if (mpn_cmp(rp1, rp2, m + n) != 0 || ret1 != ret2)
            TEST_FUNCTION_FAIL(
                    "ix = %wd\n"
                    "(m, n) = (%wd, %wd)\n"
                    "xp = %{ulong*}\n"
                    "yp = %{ulong*}\n"
                    "Expected ret: %{ulong}\n"
                    "Got ret:      %{ulong}\n"
                    "Expected rp: %{ulong*}\n"
                    "Got rp:      %{ulong*}\n",
                    ix, m, n, xp, m, yp, n, ret1, ret2, rp1, m + n, rp2, m + n);

        flint_free(xp);
        flint_free(yp);
        flint_free(rp1);
        flint_free(rp2);
    }

    TEST_FUNCTION_END(state);
}

#undef N_MIN
#undef N_MAX
#undef WANT_BIG
#undef N_MIN_STOR
#undef N_MAX_STOR
