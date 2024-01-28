/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_mod_mat.h"

TEST_FUNCTION_START(fmpz_mod_mat_scalar_mul_si, state)
{
    slong m, n, rep;

    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        fmpz_mod_mat_t A, B, C, D;
        ulong c;
        fmpz_t c1;
        fmpz_mod_ctx_t ctx;

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        fmpz_init(c1);
        fmpz_mod_ctx_init_rand_bits(ctx, state, 200);

        c = n_randtest(state);

        fmpz_mod_mat_init(A, m, n, ctx);
        fmpz_mod_mat_init(B, m, n, ctx);
        fmpz_mod_mat_init(C, m, n, ctx);
        fmpz_mod_mat_init(D, m, n, ctx);

        fmpz_mod_mat_randtest(A, state, ctx);
        fmpz_mod_mat_randtest(B, state, ctx);

        fmpz_mod_mat_scalar_mul_si(C, A, c, ctx);
        fmpz_set_si(c1, c);
        fmpz_sub_si(c1, c1, 1);
        fmpz_mod(c1, c1, ctx->n);
        fmpz_mod_mat_scalar_mul_fmpz(D, A, c1, ctx);

        /* c*A - (c-1)*A == A */
        fmpz_mod_mat_sub(D, C, D, ctx);

        if (!fmpz_mod_mat_equal(A, D, ctx))
        {
            flint_printf("FAIL\n");
            fflush(stdout);
            flint_abort();
        }

        /* Aliasing */
        fmpz_mod_mat_scalar_mul_si(C, A, c, ctx);
        fmpz_mod_mat_scalar_mul_si(A, A, c, ctx);

        if (!fmpz_mod_mat_equal(A, C, ctx))
        {
            flint_printf("FAIL\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_mat_clear(A, ctx);
        fmpz_mod_mat_clear(B, ctx);
        fmpz_mod_mat_clear(C, ctx);
        fmpz_mod_mat_clear(D, ctx);
        fmpz_mod_ctx_clear(ctx);
        fmpz_clear(c1);
    }

    TEST_FUNCTION_END(state);
}
