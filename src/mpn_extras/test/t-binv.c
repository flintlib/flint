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

TEST_FUNCTION_START(flint_mpn_binv, state)
{
    slong iter;

    for (iter = 0; iter < 5000 * flint_test_multiplier(); iter++)
    {
        mp_ptr x, y, c;
        mp_size_t xn, n;

        xn = 1 + n_randint(state, 40);
        n = 1 + n_randint(state, 60);
        if (n_randint(state, 50) == 0)
        {
            xn = 1 + n_randint(state, 1500);
            n = 1 + n_randint(state, 1500);
        }

        x = flint_malloc(xn * sizeof(mp_limb_t));
        y = flint_malloc(n * sizeof(mp_limb_t));
        c = flint_malloc(n * sizeof(mp_limb_t));

        flint_mpn_rrandom(x, state, xn);
        x[0] |= 1;

        flint_mpn_binv(y, x, xn, n);
        flint_mpn_mulmid(c, x, FLINT_MIN(xn, n), y, n, 0, n);

        if (c[0] != 1 || !flint_mpn_zero_p(c + 1, n - 1))
            TEST_FUNCTION_FAIL("xn = %wd, n = %wd\nx = %{ulong*}\ny = %{ulong*}\nxy = %{ulong*}\n",
                xn, n, x, xn, y, n, c, n);

        flint_free(x);
        flint_free(y);
        flint_free(c);
    }

    TEST_FUNCTION_END(state);
}
