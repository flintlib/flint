/*
    Copyright (C) 2025 Lars Göttgens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mod_mat.h"

TEST_FUNCTION_START(fmpz_mod_mat_transpose, state)
{
    slong m, n, rep;

    /* Rectangular transpose */
    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        fmpz_mod_mat_t A, B, C;
        fmpz_mod_ctx_t ctx;

        fmpz_mod_ctx_init_rand_bits(ctx, state, 200);

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        fmpz_mod_mat_init(A, m, n, ctx);
        fmpz_mod_mat_init(B, n, m, ctx);
        fmpz_mod_mat_init(C, m, n, ctx);

        fmpz_mod_mat_randtest(A, state, ctx);
        fmpz_mod_mat_randtest(B, state, ctx);

        fmpz_mod_mat_transpose(B, A, ctx);
        fmpz_mod_mat_transpose(C, B, ctx);

        if (!fmpz_mod_mat_equal(C, A, ctx))
        {
            flint_printf("FAIL: C != A\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_mat_clear(A, ctx);
        fmpz_mod_mat_clear(B, ctx);
        fmpz_mod_mat_clear(C, ctx);
        fmpz_mod_ctx_clear(ctx);
    }

    /* Self-transpose */
    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        fmpz_mod_mat_t A, B;
        fmpz_mod_ctx_t ctx;

        fmpz_mod_ctx_init_rand_bits(ctx, state, 200);
        
        m = n_randint(state, 20);

        fmpz_mod_mat_init(A, m, m, ctx);
        fmpz_mod_mat_init(B, m, m, ctx);

        fmpz_mod_mat_randtest(A, state, ctx);
        fmpz_mod_mat_set(B, A, ctx);
        fmpz_mod_mat_transpose(B, B, ctx);
        fmpz_mod_mat_transpose(B, B, ctx);

        if (!fmpz_mod_mat_equal(B, A, ctx))
        {
            flint_printf("FAIL: B != A\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_mat_clear(A, ctx);
        fmpz_mod_mat_clear(B, ctx);
        fmpz_mod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
