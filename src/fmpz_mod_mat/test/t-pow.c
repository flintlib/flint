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

TEST_FUNCTION_START(fmpz_mod_mat_pow, state)
{
    slong i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mat_t A, B, C;
        slong n, i;
        ulong e;
        fmpz_mod_ctx_t ctx;

        n = n_randint(state, 10);
        e = n_randint(state, 20);

        fmpz_mod_ctx_init_rand_bits(ctx, state, 200);

        fmpz_mod_mat_init(A, n, n, ctx);
        fmpz_mod_mat_init(B, n, n, ctx);
        fmpz_mod_mat_init(C, n, n, ctx);

        fmpz_mod_mat_randtest(A, state, ctx);
        fmpz_mod_mat_randtest(B, state, ctx);

        /* Make sure noise in the output is ok */
        fmpz_mod_mat_randtest(B, state, ctx);

        fmpz_mod_mat_pow_ui(B, A, e, ctx);

        fmpz_mod_mat_one(C, ctx);
        for (i = 0; i < e; i++)
            fmpz_mod_mat_mul(C, C, A, ctx);

        if (!fmpz_mod_mat_equal(C, B, ctx))
        {
            flint_printf("FAIL: results not equal\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_mat_pow_ui(A, A, e, ctx);

        if (!fmpz_mod_mat_equal(A, B, ctx))
        {
            flint_printf("FAIL: aliasing failed\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_mat_clear(A, ctx);
        fmpz_mod_mat_clear(B, ctx);
        fmpz_mod_mat_clear(C, ctx);
        fmpz_mod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
