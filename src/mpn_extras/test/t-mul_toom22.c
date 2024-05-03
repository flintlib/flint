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

/* We want m > n / 2 + 2 and m <= n. This is only possible if n > 4. */
#define N_MIN 5
#define N_MAX 300
#define M_MIN(n) ((n) / 2 + 3)
#define M_MAX(n) (n)

TEST_FUNCTION_START(flint_mpn_mul_toom22, state)
{
    slong iter;

    for (iter = 0; iter < 5000 * flint_test_multiplier(); iter++)
    {
        slong n, m;
        mp_ptr X, Y, R1, R2;

        n = N_MIN + n_randint(state, N_MAX - N_MIN + 1);
        m = M_MIN(n) + n_randint(state, M_MAX(n) - M_MIN(n) + 1);

        X = flint_malloc(sizeof(mp_limb_t) * n);
        Y = flint_malloc(sizeof(mp_limb_t) * m);
        R1 = flint_malloc(sizeof(mp_limb_t) * (n + m));
        R2 = flint_malloc(sizeof(mp_limb_t) * (n + m));

        flint_mpn_rrandom(X, state, n);
        flint_mpn_rrandom(Y, state, m);
        flint_mpn_rrandom(R1, state, n + m);

        flint_mpn_mul_toom22(R1, X, n, Y, m, NULL);
        mpn_mul(R2, X, n, Y, m);

        if (mpn_cmp(R1, R2, n + m) != 0)
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

#undef N_MIN
#undef N_MAX
#undef M_MIN
#undef M_MAX
