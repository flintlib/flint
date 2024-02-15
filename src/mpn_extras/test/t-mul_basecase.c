/*
    Copyright (C) 2023, 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"

#if FLINT_HAVE_ADX

TEST_FUNCTION_START(flint_mpn_mul_basecase, state)
{
    slong ix;

    for (ix = 0; ix < 100000 * flint_test_multiplier(); ix++)
    {
        mp_limb_t res1[32] = {UWORD(0)};
        mp_limb_t res2[32] = {UWORD(0)};
        mp_limb_t ret1, ret2;
        mp_limb_t ap[16];
        mp_limb_t bp[16];
        slong alen, blen;

        alen = 1 + n_randint(state, 16);
        blen = 1 + n_randint(state, FLINT_MIN(alen, WORD(16)));

        mpn_random2(ap, alen);
        mpn_random2(bp, blen);

        ret1 = flint_mpn_mul_basecase(res1, ap, bp, alen, blen);
        ret2 = mpn_mul(res2, ap, alen, bp, blen);

        if (mpn_cmp(res1, res2, alen + blen) || ret1 != ret2)
            TEST_FUNCTION_FAIL(
                    "ix = %wd\n"
                    "alen = %wd\n"
                    "blen = %wd\n"
                    "ap = %{ulong*}\n"
                    "bp = %{ulong*}\n"
                    "ret1 = %wu\n"
                    "ret2 = %wu\n"
                    "Got:      %{ulong*}\n"
                    "Expected: %{ulong*}\n",
                    ix, alen, blen, ap, alen, bp, blen, ret1, ret2, res1, alen + blen, res2, alen + blen);
    }

    TEST_FUNCTION_END(state);
}
#else
TEST_FUNCTION_START(flint_mpn_mul_basecase, state)
{
    TEST_FUNCTION_END_SKIPPED(state);
}
#endif
