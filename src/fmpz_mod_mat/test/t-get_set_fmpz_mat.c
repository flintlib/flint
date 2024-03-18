/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mod_mat.h"

TEST_FUNCTION_START(fmpz_mod_mat_get_set_fmpz_mat, state)
{
    slong m, n, rep;

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        fmpz_mod_mat_t A;
        fmpz_mod_mat_t B;
        fmpz_mod_ctx_t ctx;
        fmpz_mat_t C;

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        fmpz_mod_ctx_init_rand_bits(ctx, state, 200);

        fmpz_mod_mat_init(A, m, n, ctx);
        fmpz_mod_mat_init(B, m, n, ctx);

        fmpz_mat_init(C, m, n);

        fmpz_mod_mat_randtest(A, state, ctx);

        fmpz_mod_mat_get_fmpz_mat(C, A, ctx);
        fmpz_mod_mat_set_fmpz_mat(B, C, ctx);

        if (!fmpz_mod_mat_equal(A, B, ctx))
        {
            flint_printf("FAIL: matrices not equal!\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(C);

        fmpz_mod_mat_clear(A, ctx);
        fmpz_mod_mat_clear(B, ctx);
        fmpz_mod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
