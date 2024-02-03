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

TEST_FUNCTION_START(flint_mpn_mul, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * flint_test_multiplier(); iter++)
    {
        slong i, n, m;
        mp_ptr X, Y, R1, R2;
        mp_limb_t ret1, ret2;

        n = 1 + n_randint(state, 15);
        if (n_randint(state, 10000) == 0)
            n = 1 + n_randint(state, 1000);
        m = 1 + n_randint(state, n);

        X = flint_malloc(sizeof(mp_limb_t) * n);
        Y = flint_malloc(sizeof(mp_limb_t) * m);
        R1 = flint_malloc(sizeof(mp_limb_t) * (n + m));
        R2 = flint_malloc(sizeof(mp_limb_t) * (n + m));

        mpn_random2(X, n);
        mpn_random2(Y, m);

        for (i = 0; i < n + m; i++)
            R1[i] = n_randtest(state);

        ret1 = flint_mpn_mul(R1, X, n, Y, m);
        ret2 = mpn_mul(R2, X, n, Y, m);

        if (mpn_cmp(R1, R2, n + m) != 0 || ret1 != ret2)
            TEST_FUNCTION_FAIL(
                    "n = %wd\n"
                    "X = %{ulong*}\n"
                    "Y = %{ulong*}\n"
                    "R1 = %{ulong*}\n"
                    "R2 = %{ulong*}\n",
                    n, X, n, Y, m, R1, n + m, R2, n + m);

        flint_free(X);
        flint_free(Y);
        flint_free(R1);
        flint_free(R2);
    }

    TEST_FUNCTION_END(state);
}
