/*
    Copyright (C) 2010 Fredrik Johansson
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

TEST_FUNCTION_START(fmpz_mod_mat_rank, state)
{
    fmpz_mod_ctx_t ctx;
    fmpz_mod_mat_t A;
    slong i, m, n, d, r;

    /* Maximally sparse matrices of given rank */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 20);
        n = n_randint(state, 20);
        fmpz_mod_ctx_init_rand_bits_prime(ctx, state, 200);

        for (r = 0; r <= FLINT_MIN(m, n); r++)
        {
            fmpz_mod_mat_init(A, m, n, fmpz_mod_ctx_modulus(ctx));
            fmpz_mod_mat_randrank(A, state, r);
            FLINT_TEST(r == fmpz_mod_mat_rank(A));
            fmpz_mod_mat_clear(A);
        }

        fmpz_mod_ctx_clear(ctx);
    }

    /* Dense */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 20);
        n = n_randint(state, 20);
        fmpz_mod_ctx_init_rand_bits_prime(ctx, state, 200);

        for (r = 0; r <= FLINT_MIN(m, n); r++)
        {
            d = n_randint(state, 2 * m * n + 1);
            fmpz_mod_mat_init(A, m, n, fmpz_mod_ctx_modulus(ctx));
            fmpz_mod_mat_randrank(A, state, r);
            fmpz_mod_mat_randops(A, d, state);
            FLINT_TEST(r == fmpz_mod_mat_rank(A));
            fmpz_mod_mat_clear(A);
        }
        fmpz_mod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
