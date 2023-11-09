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

        fmpz_mod_mat_init(A, m, k, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_init(B, m, n, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_init(X, k, n, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_init(AX, m, n, fmpz_mod_ctx_modulus(ctx));

        fmpz_mod_mat_randrank(A, state, n_randint(state, FLINT_MIN(m, k) + 1));
        fmpz_mod_mat_randtest(B, state);

        /* Dense */
        if (n_randint(state, 2))
            fmpz_mod_mat_randops(A, 1+n_randint(state, 1+m*m), state);

        solved = fmpz_mod_mat_can_solve(X, A, B);
        fmpz_mod_mat_mul(AX, A, X);
        FLINT_TEST(!solved || fmpz_mod_mat_equal(AX, B));

        fmpz_mod_mat_clear(A);
        fmpz_mod_mat_clear(B);
        fmpz_mod_mat_clear(X);
        fmpz_mod_mat_clear(AX);

        fmpz_mod_ctx_clear(ctx);
    }

    /* test random solvable systems */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mod_ctx_init_rand_bits_prime(ctx, state, 200);

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        k = n_randint(state, 20);

        fmpz_mod_mat_init(A, m, k, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_init(B, m, n, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_init(X, k, n, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_init(X2, k, n, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_init(AX, m, n, fmpz_mod_ctx_modulus(ctx));

        fmpz_mod_mat_randrank(A, state, n_randint(state, FLINT_MIN(m, k) + 1));
        fmpz_mod_mat_randtest(X2, state);

        fmpz_mod_mat_mul(B, A, X2);
        solved = fmpz_mod_mat_can_solve(X, A, B);
        fmpz_mod_mat_mul(AX, A, X);
        FLINT_TEST(solved && fmpz_mod_mat_equal(AX, B));

        fmpz_mod_mat_clear(A);
        fmpz_mod_mat_clear(B);
        fmpz_mod_mat_clear(X);
        fmpz_mod_mat_clear(X2);
        fmpz_mod_mat_clear(AX);

        fmpz_mod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
