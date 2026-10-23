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

TEST_FUNCTION_START(flint_mpn_brsqrt, state)
{
    slong iter;

    for (iter = 0; iter < 3000 * flint_test_multiplier(); iter++)
    {
        mp_ptr x, y, c, d;
        mp_size_t xn, n;
        int res;

        xn = 1 + n_randint(state, 40);
        n = 1 + n_randint(state, 60);
        if (n_randint(state, 50) == 0)
        {
            xn = 1 + n_randint(state, 1200);
            n = 1 + n_randint(state, 1200);
        }

        x = flint_malloc(xn * sizeof(mp_limb_t));
        y = flint_malloc(n * sizeof(mp_limb_t));
        c = flint_malloc(n * sizeof(mp_limb_t));
        d = flint_malloc(n * sizeof(mp_limb_t));

        flint_mpn_rrandom(x, state, xn);
        if (n_randint(state, 4) != 0)
            x[0] = (x[0] & ~UWORD(7)) | 1;

        res = flint_mpn_brsqrt(y, x, xn, n);

        if (res != ((x[0] & 7) == 1))
            TEST_FUNCTION_FAIL("wrong return value: xn = %wd, n = %wd, x0 = %wu\n", xn, n, x[0]);

        if (res)
        {
            /* x y^2 == 1 mod B^n, y == 1 mod 4 */
            flint_mpn_mulmid(c, y, n, y, n, 0, n);
            flint_mpn_mulmid(d, c, n, x, FLINT_MIN(xn, n), 0, n);

            if (d[0] != 1 || !flint_mpn_zero_p(d + 1, n - 1) || (y[0] & 3) != 1)
                TEST_FUNCTION_FAIL("xn = %wd, n = %wd\nx = %{ulong*}\ny = %{ulong*}\nx y^2 = %{ulong*}\n",
                    xn, n, x, xn, y, n, d, n);
        }

        flint_free(x);
        flint_free(y);
        flint_free(c);
        flint_free(d);
    }

    TEST_FUNCTION_END(state);
}
