/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
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

        fmpz_mod_mat_init(A, m, m, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_init(B, m, m, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_init(C, m, m, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_init(I, m, m, fmpz_mod_ctx_modulus(ctx));

        for (j = 0; j < m; j++)
            fmpz_one(fmpz_mod_mat_entry(I, j, j));

        /* Verify that A * A^-1 = I for random matrices */

        fmpz_mod_mat_randrank(A, state, m);
        /* Dense or sparse? */
        if (n_randint(state, 2))
            fmpz_mod_mat_randops(A, 1+n_randint(state, 1+m*m), state);

        result = fmpz_mod_mat_inv(B, A);
        fmpz_mod_mat_mul(C, A, B);
        FLINT_TEST(result && fmpz_mod_mat_equal(C, I));

        /* Test aliasing */
        fmpz_mod_mat_set(C, A);
        fmpz_mod_mat_inv(A, A);
        fmpz_mod_mat_mul(B, A, C);
        FLINT_TEST(fmpz_mod_mat_equal(B, I));

        fmpz_mod_mat_clear(A);
        fmpz_mod_mat_clear(B);
        fmpz_mod_mat_clear(C);
        fmpz_mod_mat_clear(I);

        fmpz_mod_ctx_clear(ctx);
    }

    /* Test singular systems */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = 1 + n_randint(state, 20);
        r = n_randint(state, m);
        fmpz_mod_ctx_init_rand_bits_prime(ctx, state, 200);

        fmpz_mod_mat_init(A, m, m, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_mat_init(B, m, m, fmpz_mod_ctx_modulus(ctx));

        fmpz_mod_mat_randrank(A, state, r);

        /* Dense */
        if (n_randint(state, 2))
            fmpz_mod_mat_randops(A, 1+n_randint(state, 1+m*m), state);

        FLINT_TEST(!fmpz_mod_mat_inv(B, A));
        FLINT_TEST(!fmpz_mod_mat_inv(A, A));

        fmpz_mod_mat_clear(A);
        fmpz_mod_mat_clear(B);

        fmpz_mod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
