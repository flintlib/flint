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

TEST_FUNCTION_START(flint_mpn_inv, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        mp_ptr x, q, q2, N, r;
        mp_size_t xn, n, qn;

        if (n_randint(state, 3) == 0)
        {
            xn = 1 + n_randint(state, 1000);
            n = xn + n_randint(state, 1200);
        }
        else
        {
            xn = 1 + n_randint(state, 20);
            n = xn + n_randint(state, 30);
        }
        qn = n - xn + 2;

        x = flint_malloc(xn * sizeof(mp_limb_t));
        q = flint_malloc(qn * sizeof(mp_limb_t));
        q2 = flint_malloc(qn * sizeof(mp_limb_t));
        N = flint_malloc((n + 1) * sizeof(mp_limb_t));
        r = flint_malloc(xn * sizeof(mp_limb_t));

        flint_mpn_rrandom(x, state, xn);
        if (x[xn - 1] == 0)
            x[xn - 1] = 1;
        if (n_randint(state, 20) == 0)
        {
            flint_mpn_zero(x, xn);
            x[xn - 1] = 1;
        }

        flint_mpn_inv(q, x, xn, n);

        flint_mpn_zero(N, n);
        N[n] = 1;
        mpn_tdiv_qr(q2, r, 0, N, n + 1, x, xn);

        if (mpn_cmp(q, q2, qn) != 0)
            TEST_FUNCTION_FAIL("xn = %wd, n = %wd\nx = %{ulong*}\nq = %{ulong*}\nq2 = %{ulong*}\n",
                xn, n, x, xn, q, qn, q2, qn);

        flint_free(x);
        flint_free(q);
        flint_free(q2);
        flint_free(N);
        flint_free(r);
    }

    TEST_FUNCTION_END(state);
}
