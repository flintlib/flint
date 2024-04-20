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

#if FLINT_HAVE_NATIVE_mpn_2add_n_inplace

#define N_MIN 4
#define N_MAX 100

TEST_FUNCTION_START(flint_mpn_2add_n_inplace, state)
{
    slong ix;
    int result;

    for (ix = 0; ix < 100000 * flint_test_multiplier(); ix++)
    {
        mp_limb_t rp[N_MAX];
        mp_limb_t ap[N_MAX];
        mp_limb_t bp[N_MAX];
        mp_limb_t rp1[N_MAX];
        mp_limb_t rp2[N_MAX];
        mp_limb_t c1, c2;
        mp_size_t n;

        n = N_MIN + n_randint(state, N_MAX - N_MIN + 1);

        flint_mpn_rrandom(rp, state, N_MAX);
        flint_mpn_rrandom(ap, state, N_MAX);
        flint_mpn_rrandom(bp, state, N_MAX);
        flint_mpn_copyi(rp1, rp, N_MAX);
        flint_mpn_copyi(rp2, rp, N_MAX);

        c2 = flint_mpn_2add_n_inplace(rp2, ap, bp, n);

        c1 = mpn_add_n(rp1, rp1, ap, n);
        c1 += mpn_add_n(rp1, rp1, bp, n);

        result = (c1 == c2 && mpn_cmp(rp1, rp2, N_MAX) == 0);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "ix = %wd\n"
                    "n = %wd\n"
                    "ap = %{ulong*}\n"
                    "bp = %{ulong*}\n"
                    "rp = %{ulong*}\n"
                    "Expected carry: %{ulong}\n"
                    "Got carry:      %{ulong}\n"
                    "Expected limbs: %{ulong*}\n"
                    "Got limbs:      %{ulong*}\n",
                    ix, n, ap, n, bp, n, rp, n, c1, c2, rp1, n, rp2, n);
    }

    TEST_FUNCTION_END(state);
}
# undef N_MIN
# undef N_MAX
#else
TEST_FUNCTION_START(flint_mpn_2add_n_inplace, state)
{
    TEST_FUNCTION_END_SKIPPED(state);
}
#endif
