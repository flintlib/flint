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

#if FLINT_HAVE_ADX

TEST_FUNCTION_START(flint_mpn_sqr_basecase, state)
{
    slong ix;

    for (ix = 0; ix < 100000 * flint_test_multiplier(); ix++)
    {
        mp_limb_t res1[14] = {UWORD(0)};
        mp_limb_t res2[14] = {UWORD(0)};
        mp_limb_t ret1;
        mp_limb_t ap[7];
        slong alen;

        alen = 1 + n_randint(state, 7);

        mpn_random2(ap, alen);

        ret1 = flint_mpn_sqr_basecase(res1, ap, alen);
        mpn_sqr(res2, ap, alen);

        if (mpn_cmp(res1, res2, 2 * alen) || ret1 != res1[2 * alen - 1])
            TEST_FUNCTION_FAIL(
                    "ix = %wd\n"
                    "alen = %wd\n"
                    "ap = %{ulong*}\n"
                    "ret1 = %wu\n"
                    "Got:      %{ulong*}\n"
                    "Expected: %{ulong*}\n",
                    ix, alen, ap, alen, ret1, res1, 2 * alen, res2, 2 * alen);
    }

    TEST_FUNCTION_END(state);
}
#else
TEST_FUNCTION_START(flint_mpn_sqr_basecase, state)
{
    TEST_FUNCTION_END_SKIPPED(state);
}
#endif
