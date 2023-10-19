/*
    Copyright (C) 2010, 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2018 Tommy Hofmann
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

TEST_FUNCTION_START(fmpz_mod_mat_solve, state)
{
    fmpz_mod_ctx_t ctx;
    fmpz_mod_mat_t A, X, B, AX;
    slong i, m, n, r;
    int solved;

    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        fmpz_mod_ctx_init_rand_bits_prime(ctx, state, 200);

        m = n_randint(state, 50);
        n = n_randint(state, 50);

        fmpz_mod_mat_init(A, m, m, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_init(B, m, n, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_init(X, m, n, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_init(AX, m, n, fmpz_mod_ctx_modulus(ctx));

        fmpz_mod_mat_randrank(A, state, m);
        fmpz_mod_mat_randtest(B, state);

        /* Dense */
        if (n_randint(state, 2))
            fmpz_mod_mat_randops(A, 1+n_randint(state, 1+m*m), state);

        solved = fmpz_mod_mat_solve(X, A, B);
        fmpz_mod_mat_mul(AX, A, X);
        FLINT_TEST(solved && fmpz_mod_mat_equal(AX, B));

        fmpz_mod_mat_clear(A);
        fmpz_mod_mat_clear(B);
        fmpz_mod_mat_clear(X);
        fmpz_mod_mat_clear(AX);

        fmpz_mod_ctx_clear(ctx);
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mod_ctx_init_rand_bits_prime(ctx, state, 200);

        m = 1 + n_randint(state, 20);
        n = 1 + n_randint(state, 20);
        r = n_randint(state, m);

        fmpz_mod_mat_init(A, m, m, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_init(B, m, n, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_init(X, m, n, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_init(AX, m, n, fmpz_mod_ctx_modulus(ctx));

        fmpz_mod_mat_randrank(A, state, r);
        fmpz_mod_mat_randtest(B, state);

        /* Dense */
        if (n_randint(state, 2))
            fmpz_mod_mat_randops(A, 1+n_randint(state, 1+m*m), state);

        solved = fmpz_mod_mat_solve(X, A, B);
        FLINT_TEST(!solved);

        fmpz_mod_mat_clear(A);
        fmpz_mod_mat_clear(B);
        fmpz_mod_mat_clear(X);
        fmpz_mod_mat_clear(AX);

        fmpz_mod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
