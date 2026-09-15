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

TEST_FUNCTION_START(fmpz_divmod_2exp, state)
{
    slong i;

    for (i = 0; i < 5000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, q, c, m, y;
        flint_bitcnt_t N;

        fmpz_init(a); fmpz_init(b); fmpz_init(q); fmpz_init(c); fmpz_init(m); fmpz_init(y);

        N = 1 + n_randint(state, n_randint(state, 20) == 0 ? 100000 : 400);
        fmpz_randtest(a, state, 1 + n_randint(state, N + 200));
        fmpz_randtest(b, state, 1 + n_randint(state, N + 200));
        fmpz_setbit(b, 0);
        fmpz_one(m);
        fmpz_mul_2exp(m, m, N);

        fmpz_divmod_2exp(q, a, b, N);
        fmpz_mul(c, q, b);
        fmpz_sub(c, c, a);
        fmpz_mod(c, c, m);

        if (!fmpz_is_zero(c) || fmpz_sgn(q) < 0 || fmpz_cmp(q, m) >= 0)
            TEST_FUNCTION_FAIL("N = %wd\na = %{fmpz}\nb = %{fmpz}\nq = %{fmpz}\n", (slong) N, a, b, q);

        /* rsqrt: y^2 a == 1 mod 2^N when a == 1 mod 8 */
        fmpz_fdiv_q_2exp(a, a, 3);
        fmpz_mul_2exp(a, a, 3);
        fmpz_add_ui(a, a, 1);
        if (!fmpz_rsqrtmod_2exp(y, a, N))
            TEST_FUNCTION_FAIL("rsqrtmod_2exp returned 0: N = %wd\n", (slong) N);
        fmpz_mul(c, y, y);
        fmpz_mul(c, c, a);
        fmpz_sub_ui(c, c, 1);
        fmpz_mod(c, c, m);
        if (!fmpz_is_zero(c) || fmpz_sgn(y) < 0 || fmpz_cmp(y, m) >= 0)
            TEST_FUNCTION_FAIL("rsqrtmod_2exp: N = %wd\na = %{fmpz}\ny = %{fmpz}\n", (slong) N, a, y);

        fmpz_clear(a); fmpz_clear(b); fmpz_clear(q); fmpz_clear(c); fmpz_clear(m); fmpz_clear(y);
    }

    TEST_FUNCTION_END(state);
}
