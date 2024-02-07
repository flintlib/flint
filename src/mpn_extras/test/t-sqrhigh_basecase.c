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
#if FLINT_HAVE_NATIVE_MPN_SQRHIGH_BASECASE

# define N_MIN 1
# define N_MAX 64

TEST_FUNCTION_START(flint_mpn_sqrhigh_basecase, state)
{
    slong ix;
    int result;

    for (ix = 0; ix < 100000 * flint_test_multiplier(); ix++)
    {
        mp_limb_t rp1[N_MAX + 1];
        mp_limb_t rp2[N_MAX + 1];
        mp_limb_t xp[N_MAX];
        mp_size_t n;

        n = N_MIN + n_randint(state, N_MAX - N_MIN + 1);

        mpn_random2(xp, n);

        rp1[0] = flint_mpn_mulhigh_basecase(rp1 + 1, xp, xp, n);
        rp2[0] = flint_mpn_sqrhigh_basecase(rp2 + 1, xp, n);

        result = (mpn_cmp(rp1, rp2, n + 1) == 0);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "Wrong result!\n"
                    "ix = %wd\n"
                    "n = %wd\n"
                    "xp = %{ulong*}\n"
                    "Expected: %{ulong*}\n"
                    "Got:      %{ulong*}\n",
                    ix, n, xp, n, rp1, n + 1, rp2, n + 1);
    }

    TEST_FUNCTION_END(state);
}
# undef N_MIN
# undef N_MAX
#else
TEST_FUNCTION_START(flint_mpn_sqrhigh_basecase, state)
{
    TEST_FUNCTION_END_SKIPPED(state);
}
#endif
