/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mat.h"

TEST_FUNCTION_START(fmpz_mod_mat_inv, state)
{
    fmpz_mod_mat_t A, B, C, I;
    fmpz_mod_ctx_t ctx;
    slong i, j, m, r;
    int result;

    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 20);

        fmpz_mod_ctx_init_rand_bits_prime(ctx, state, 200);

        fmpz_mod_mat_init(A, m, m, ctx);
        fmpz_mod_mat_init(B, m, m, ctx);
        fmpz_mod_mat_init(C, m, m, ctx);
        fmpz_mod_mat_init(I, m, m, ctx);

        for (j = 0; j < m; j++)
            fmpz_one(fmpz_mod_mat_entry(I, j, j));

        /* Verify that A * A^-1 = I for random matrices */

        fmpz_mod_mat_randrank(A, state, m, ctx);
        /* Dense or sparse? */
        if (n_randint(state, 2))
            fmpz_mod_mat_randops(A, state, 1+n_randint(state, 1+m*m), ctx);

        result = fmpz_mod_mat_inv(B, A, ctx);
        fmpz_mod_mat_mul(C, A, B, ctx);
        FLINT_TEST(result && fmpz_mod_mat_equal(C, I, ctx));

        /* Test aliasing */
        fmpz_mod_mat_set(C, A, ctx);
        fmpz_mod_mat_inv(A, A, ctx);
        fmpz_mod_mat_mul(B, A, C, ctx);
        FLINT_TEST(fmpz_mod_mat_equal(B, I, ctx));

        fmpz_mod_mat_clear(A, ctx);
        fmpz_mod_mat_clear(B, ctx);
        fmpz_mod_mat_clear(C, ctx);
        fmpz_mod_mat_clear(I, ctx);

        fmpz_mod_ctx_clear(ctx);
    }

    /* Test singular systems */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = 1 + n_randint(state, 20);
        r = n_randint(state, m);
        fmpz_mod_ctx_init_rand_bits_prime(ctx, state, 200);

        fmpz_mod_mat_init(A, m, m, ctx);
        fmpz_mod_mat_init(B, m, m, ctx);

        fmpz_mod_mat_randrank(A, state, r, ctx);

        /* Dense */
        if (n_randint(state, 2))
            fmpz_mod_mat_randops(A, state, 1+n_randint(state, 1+m*m), ctx);

        FLINT_TEST(!fmpz_mod_mat_inv(B, A, ctx));
        FLINT_TEST(!fmpz_mod_mat_inv(A, A, ctx));

        fmpz_mod_mat_clear(A, ctx);
        fmpz_mod_mat_clear(B, ctx);

        fmpz_mod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
