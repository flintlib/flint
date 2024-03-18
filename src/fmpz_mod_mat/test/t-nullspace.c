/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mat.h"

TEST_FUNCTION_START(fmpz_mod_mat_nullspace, state)
{
    slong i;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mod_ctx_t ctx;
        fmpz_mod_mat_t A, B, ker;
        slong m, n, d, r, nullity, nulrank;

        m = n_randint(state, 30);
        n = n_randint(state, 30);

        for (r = 0; r <= FLINT_MIN(m, n); r++)
        {
            fmpz_mod_ctx_init_rand_bits_prime(ctx, state, 200);
            d = n_randint(state, 2 * m * n + 1);

            fmpz_mod_mat_init(A, m, n, ctx);
            fmpz_mod_mat_init(ker, n, n, ctx);
            fmpz_mod_mat_init(B, m, n, ctx);

            fmpz_mod_mat_randrank(A, state, r, ctx);
            /* Densify */
            if (n_randlimb(state) % 2)
                fmpz_mod_mat_randops(A, state, d, ctx);

            nullity = fmpz_mod_mat_nullspace(ker, A, ctx);
            nulrank = fmpz_mod_mat_rank(ker, ctx);
            FLINT_TEST(nullity == nulrank);
            FLINT_TEST(nullity + r == n);

            fmpz_mod_mat_mul(B, A, ker, ctx);
            FLINT_TEST(fmpz_mod_mat_rank(B, ctx) == 0);

            fmpz_mod_mat_clear(A, ctx);
            fmpz_mod_mat_clear(ker, ctx);
            fmpz_mod_mat_clear(B, ctx);

            fmpz_mod_ctx_clear(ctx);
        }
    }

    TEST_FUNCTION_END(state);
}
