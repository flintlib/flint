/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"

TEST_FUNCTION_START(flint_mpn_bsqrt, state)
{
    slong iter;

    for (iter = 0; iter < 3000 * flint_test_multiplier(); iter++)
    {
        mp_ptr x, s, c;
        mp_size_t xn, n, i;
        int res;

        xn = 1 + n_randint(state, 40);
        n = 1 + n_randint(state, 60);
        if (n_randint(state, 50) == 0)
        {
            xn = 1 + n_randint(state, 1200);
            n = 1 + n_randint(state, 1200);
        }

        x = flint_malloc(xn * sizeof(mp_limb_t));
        s = flint_malloc(n * sizeof(mp_limb_t));
        c = flint_malloc(n * sizeof(mp_limb_t));

        flint_mpn_rrandom(x, state, xn);
        if (n_randint(state, 4) != 0)
            x[0] = (x[0] & ~UWORD(7)) | 1;

        res = flint_mpn_bsqrt(s, x, xn, n);

        if (res != ((x[0] & 7) == 1))
            TEST_FUNCTION_FAIL("wrong return value: xn = %wd, n = %wd, x0 = %wu\n", xn, n, x[0]);

        if (res)
        {
            flint_mpn_mulmid(c, s, n, s, n, 0, n);

            for (i = 0; i < n; i++)
            {
                if (c[i] != ((i < xn) ? x[i] : 0))
                    TEST_FUNCTION_FAIL("xn = %wd, n = %wd\nx = %{ulong*}\ns = %{ulong*}\ns^2 = %{ulong*}\n",
                        xn, n, x, xn, s, n, c, n);
            }

            if ((s[0] & 3) != 1)
                TEST_FUNCTION_FAIL("root not 1 mod 4: xn = %wd, n = %wd\n", xn, n);
        }

        flint_free(x);
        flint_free(s);
        flint_free(c);
    }

    TEST_FUNCTION_END(state);
}
