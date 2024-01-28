/*
    Copyright (C) 2010, 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2018 Tommy Hofmann
    Copyright (C) 2020 William Hart
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

TEST_FUNCTION_START(fmpz_mod_mat_can_solve, state)
{
    fmpz_mod_ctx_t ctx;
    fmpz_mod_mat_t A, X, X2, B, AX;
    slong i, k, m, n;
    int solved;

    /* test random systems */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_mod_ctx_init_rand_bits_prime(ctx, state, 200);

        m = n_randint(state, 50);
        n = n_randint(state, 50);
        k = n_randint(state, 50);

        fmpz_mod_mat_init(A, m, k, ctx);
        fmpz_mod_mat_init(B, m, n, ctx);
        fmpz_mod_mat_init(X, k, n, ctx);
        fmpz_mod_mat_init(AX, m, n, ctx);

        fmpz_mod_mat_randrank(A, state, n_randint(state, FLINT_MIN(m, k) + 1), ctx);
        fmpz_mod_mat_randtest(B, state, ctx);

        /* Dense */
        if (n_randint(state, 2))
            fmpz_mod_mat_randops(A, state, 1+n_randint(state, 1+m*m), ctx);

        solved = fmpz_mod_mat_can_solve(X, A, B, ctx);
        fmpz_mod_mat_mul(AX, A, X, ctx);
        FLINT_TEST(!solved || fmpz_mod_mat_equal(AX, B, ctx));

        fmpz_mod_mat_clear(A, ctx);
        fmpz_mod_mat_clear(B, ctx);
        fmpz_mod_mat_clear(X, ctx);
        fmpz_mod_mat_clear(AX, ctx);

        fmpz_mod_ctx_clear(ctx);
    }

    /* test random solvable systems */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mod_ctx_init_rand_bits_prime(ctx, state, 200);

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        k = n_randint(state, 20);

        fmpz_mod_mat_init(A, m, k, ctx);
        fmpz_mod_mat_init(B, m, n, ctx);
        fmpz_mod_mat_init(X, k, n, ctx);
        fmpz_mod_mat_init(X2, k, n, ctx);
        fmpz_mod_mat_init(AX, m, n, ctx);

        fmpz_mod_mat_randrank(A, state, n_randint(state, FLINT_MIN(m, k) + 1), ctx);
        fmpz_mod_mat_randtest(X2, state, ctx);

        fmpz_mod_mat_mul(B, A, X2, ctx);
        solved = fmpz_mod_mat_can_solve(X, A, B, ctx);
        fmpz_mod_mat_mul(AX, A, X, ctx);
        FLINT_TEST(solved && fmpz_mod_mat_equal(AX, B, ctx));

        fmpz_mod_mat_clear(A, ctx);
        fmpz_mod_mat_clear(B, ctx);
        fmpz_mod_mat_clear(X, ctx);
        fmpz_mod_mat_clear(X2, ctx);
        fmpz_mod_mat_clear(AX, ctx);

        fmpz_mod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
