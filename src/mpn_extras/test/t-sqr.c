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
#define N_MAX 20

#define WANT_BIG 1
#define N_MIN_STOR (N_MAX + 1)
#define N_MAX_STOR 1000

TEST_FUNCTION_START(flint_mpn_sqr, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * flint_test_multiplier(); iter++)
    {
        slong i, n;
        mp_ptr xp, rp1, rp2;

        n = N_MIN + n_randint(state, N_MAX - N_MIN + 1);
#if WANT_BIG
        if (n_randint(state, 10000) == 0)
            n = N_MIN_STOR + n_randint(state, N_MAX_STOR - N_MIN_STOR + 1);
#endif

        xp = flint_malloc(sizeof(mp_limb_t) * n);
        rp1 = flint_malloc(sizeof(mp_limb_t) * 2 * n);
        rp2 = flint_malloc(sizeof(mp_limb_t) * 2 * n);

        flint_mpn_rrandom(xp, state, n);

        for (i = 0; i < 2 * n; i++)
            rp1[i] = n_randtest(state);

        flint_mpn_sqr(rp1, xp, n);
        mpn_sqr(rp2, xp, n);

        if (mpn_cmp(rp1, rp2, 2 * n) != 0)
            TEST_FUNCTION_FAIL(
                    "n = %wd\n"
                    "xp = %{ulong*}\n"
                    "Expected: %{ulong*}\n"
                    "Got:      %{ulong*}\n",
                    n, xp, n, rp2, 2 * n, rp1, 2 * n);

        flint_free(xp);
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
