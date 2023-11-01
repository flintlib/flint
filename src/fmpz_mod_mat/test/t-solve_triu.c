/*
    Copyright (C) 2010, 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mat.h"

TEST_FUNCTION_START(fmpz_mod_mat_solve_triu, state)
{
    slong i;

    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        fmpz_mod_ctx_t ctx;
        fmpz_mod_mat_t A, X, B, Y;
        slong rows, cols;
        int unit;

        fmpz_mod_ctx_init_rand_bits_prime(ctx, state, 200);

        rows = n_randint(state, 50);
        cols = n_randint(state, 50);
        unit = n_randint(state, 2);

        fmpz_mod_mat_init(A, rows, rows, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_init(B, rows, cols, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_init(X, rows, cols, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_init(Y, rows, cols, fmpz_mod_ctx_modulus(ctx));

        fmpz_mod_mat_randtriu(A, state, unit);
        fmpz_mod_mat_randtest(X, state);
        fmpz_mod_mat_mul(B, A, X);

        /* Check Y = A^(-1) * (A * X) = X */
        fmpz_mod_mat_solve_triu(Y, A, B, unit);
        FLINT_TEST(fmpz_mod_mat_equal(Y, X));

        /* Check aliasing */
        fmpz_mod_mat_solve_triu(B, A, B, unit);
        FLINT_TEST(fmpz_mod_mat_equal(B, X));

        fmpz_mod_mat_clear(A);
        fmpz_mod_mat_clear(B);
        fmpz_mod_mat_clear(X);
        fmpz_mod_mat_clear(Y);

        fmpz_mod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
