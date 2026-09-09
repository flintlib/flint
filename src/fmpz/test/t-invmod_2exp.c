/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_invmod_2exp, state)
{
    slong i;

    for (i = 0; i < 5000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c, m;
        flint_bitcnt_t N;

        fmpz_init(a); fmpz_init(b); fmpz_init(c); fmpz_init(m);

        N = 1 + n_randint(state, n_randint(state, 20) == 0 ? 100000 : 400);
        fmpz_randtest(a, state, 1 + n_randint(state, N + 200));
        fmpz_setbit(a, 0);

        fmpz_invmod_2exp(b, a, N);
        fmpz_one(m);
        fmpz_mul_2exp(m, m, N);
        fmpz_mul(c, a, b);
        fmpz_mod(c, c, m);

        if (!fmpz_is_one(c) || fmpz_sgn(b) < 0 || fmpz_cmp(b, m) >= 0)
            TEST_FUNCTION_FAIL("N = %wd\na = %{fmpz}\nb = %{fmpz}\n", (slong) N, a, b);

        /* aliasing */
        fmpz_set(c, a);
        fmpz_invmod_2exp(c, c, N);
        if (!fmpz_equal(b, c))
            TEST_FUNCTION_FAIL("aliasing, N = %wd\n", (slong) N);

        fmpz_clear(a); fmpz_clear(b); fmpz_clear(c); fmpz_clear(m);
    }

    TEST_FUNCTION_END(state);
}
