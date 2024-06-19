/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"

TEST_FUNCTION_START(flint_mpn_mul_n, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * flint_test_multiplier(); iter++)
    {
        slong i, n;
        mp_ptr X, Y, R1, R2;

        n = 1 + n_randint(state, 15);
        if (n_randint(state, 10000) == 0)
            n = 1 + n_randint(state, 1000);

        X = flint_malloc(sizeof(mp_limb_t) * n);
        Y = flint_malloc(sizeof(mp_limb_t) * n);
        R1 = flint_malloc(sizeof(mp_limb_t) * 2 * n);
        R2 = flint_malloc(sizeof(mp_limb_t) * 2 * n);

        flint_mpn_rrandom(X, state, n);
        flint_mpn_rrandom(Y, state, n);

        for (i = 0; i < 2 * n; i++)
            R1[i] = n_randtest(state);

        flint_mpn_mul_n(R1, X, Y, n);
        mpn_mul_n(R2, X, Y, n);

        if (mpn_cmp(R1, R2, 2 * n) != 0)
            TEST_FUNCTION_FAIL(
                    "n = %wd\n"
                    "X = %{ulong*}\n"
                    "Y = %{ulong*}\n"
                    "R1 = %{ulong*}\n"
                    "R2 = %{ulong*}\n",
                    n, X, n, Y, n, R1, 2 * n, R2, 2 * n);

        flint_free(X);
        flint_free(Y);
        flint_free(R1);
        flint_free(R2);
    }

    TEST_FUNCTION_END(state);
}
