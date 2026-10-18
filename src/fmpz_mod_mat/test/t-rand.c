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
#include "fmpz_mod_mat.h"

TEST_FUNCTION_START(fmpz_mod_mat_rand, state)
{
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        fmpz_mod_ctx_t ctx;
        fmpz_mod_mat_t A, B;
        slong r, c, i, j, k, last_zero;

        r = n_randint(state, 6);
        c = n_randint(state, 6);

        if (n_randint(state, 4) == 0)
            fmpz_mod_ctx_init_ui(ctx, 1 + n_randint(state, 3));
        else
            fmpz_mod_ctx_init_rand_bits(ctx, state, 100);

        fmpz_mod_mat_init(B, r + 1, c, ctx);

        /* all entries are reduced; the last column is not always zero */
        last_zero = 0;
        for (k = 0; k < 10; k++)
        {
            fmpz_mod_mat_zero(B, ctx);
            fmpz_mod_mat_window_init(A, B, 0, 0, r, c, ctx);
            fmpz_mod_mat_rand(A, state, ctx);
            fmpz_mod_mat_window_clear(A, ctx);

            for (i = 0; i < r; i++)
                for (j = 0; j < c; j++)
                    if (fmpz_sgn(fmpz_mod_mat_entry(B, i, j)) < 0 ||
                        fmpz_cmp(fmpz_mod_mat_entry(B, i, j), ctx->n) >= 0)
                        TEST_FUNCTION_FAIL("rand: entry out of range\n");

            for (j = 0; j < c; j++)
                if (!fmpz_is_zero(fmpz_mod_mat_entry(B, r, j)))
                    TEST_FUNCTION_FAIL("rand: wrote past the matrix\n");

            if (r > 0 && c > 0)
                last_zero += fmpz_is_zero(fmpz_mod_mat_entry(B, r - 1, c - 1));
        }

        if (r > 0 && c > 0 && fmpz_bits(ctx->n) > 20 && last_zero == 10)
            TEST_FUNCTION_FAIL("rand: last entry never set\n");

        /* randfull: entries in [1, n), or 0 if n = 1 */
        fmpz_mod_mat_zero(B, ctx);
        fmpz_mod_mat_window_init(A, B, 0, 0, r, c, ctx);
        fmpz_mod_mat_randfull(A, state, ctx);
        fmpz_mod_mat_window_clear(A, ctx);

        for (i = 0; i < r; i++)
            for (j = 0; j < c; j++)
                if (fmpz_sgn(fmpz_mod_mat_entry(B, i, j)) < 0 ||
                    fmpz_cmp(fmpz_mod_mat_entry(B, i, j), ctx->n) >= 0 ||
                    (fmpz_is_zero(fmpz_mod_mat_entry(B, i, j)) != fmpz_is_one(ctx->n)))
                    TEST_FUNCTION_FAIL("randfull: bad entry\n");

        for (j = 0; j < c; j++)
            if (!fmpz_is_zero(fmpz_mod_mat_entry(B, r, j)))
                TEST_FUNCTION_FAIL("randfull: wrote past the matrix\n");

        fmpz_mod_mat_clear(B, ctx);
        fmpz_mod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
