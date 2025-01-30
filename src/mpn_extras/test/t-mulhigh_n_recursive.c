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

TEST_FUNCTION_START(flint_mpn_mulhigh_n_recursive, state)
{
    slong ix;
    int result;
    mp_ptr rp, rpc, xp, yp;

    for (ix = 0; ix < 10000 * flint_test_multiplier(); ix++)
    {
        mp_size_t n;
        mp_limb_t ret0, ret1;

        n = 1 + n_randint(state, 100);

        xp = flint_malloc(sizeof(mp_limb_t) * n);
        yp = flint_malloc(sizeof(mp_limb_t) * n);
        rp = flint_malloc(sizeof(mp_limb_t) * n);
        rpc = flint_malloc(sizeof(mp_limb_t) * n);

        flint_mpn_rrandom(xp, state, n);
        flint_mpn_rrandom(yp, state, n);

        ret0 = _flint_mpn_mulhigh_n_naive(rpc, xp, yp, n);
        ret1 = _flint_mpn_mulhigh_n_recursive(rp, xp, yp, n);

        result = (mpn_cmp(rp, rpc, n) == 0 && ret0 == ret1);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "Not coherent with generic version.\n"
                    "ix = %wd\n"
                    "n = %wd\n"
                    "xp = %{ulong*}\n"
                    "yp = %{ulong*}\n"
                    "Expected ret: %{ulong}\n"
                    "Got ret:      %{ulong}\n"
                    "Expected limbs: %{ulong*}\n"
                    "Got limbs:      %{ulong*}\n",
                    ix, n, xp, n, yp, n, ret0, ret1, rpc, n, rp, n);

        flint_free(xp);
        flint_free(yp);
        flint_free(rp);
        flint_free(rpc);
    }

    TEST_FUNCTION_END(state);
}
