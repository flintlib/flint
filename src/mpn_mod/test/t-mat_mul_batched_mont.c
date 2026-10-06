/*
    Copyright (C) 2026 Brian Heckel

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_mod.h"
#include "gr_mat.h"
#include "fmpz.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(mpn_mod_mat_mul_batched_mont, state)
{
    slong iter;

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        fmpz_t N, val;
        gr_mat_t A, B, C, D;
        slong i, j, bits;

        /* odd modulus with 65 <= bits <= 246 (>= 2 limbs; lazy REDC valid) */
        bits = 65 + n_randint(state, 246 - 65 + 1);
        fmpz_init(N);
        fmpz_randtest_unsigned(N, state, bits);
        fmpz_setbit(N, bits - 1);
        fmpz_setbit(N, 0);

        if (gr_ctx_init_mpn_mod(ctx, N) != GR_SUCCESS)
        {
            fmpz_clear(N);
            continue;
        }

        gr_mat_init(A, 4, 4, ctx);
        gr_mat_init(B, 4, 4, ctx);
        gr_mat_init(C, 4, 4, ctx);
        gr_mat_init(D, 4, 4, ctx);

        fmpz_init(val);
        for (i = 0; i < 4; i++)
            for (j = 0; j < 4; j++)
            {
                fmpz_randm(val, state, N);
                GR_MUST_SUCCEED(gr_set_fmpz(gr_mat_entry_ptr(A, i, j, ctx), val, ctx));
                fmpz_randm(val, state, N);
                GR_MUST_SUCCEED(gr_set_fmpz(gr_mat_entry_ptr(B, i, j, ctx), val, ctx));
            }
        fmpz_clear(val);

        GR_MUST_SUCCEED(gr_mat_mul_classical(D, A, B, ctx));

        GR_MUST_SUCCEED(mpn_mod_mat_mul_batched_mont(C, A, B, ctx));
        if (gr_mat_equal(C, D, ctx) != T_TRUE)
        {
            flint_printf("FAIL (A*B): iter %wd\n", iter);
            flint_printf("n = "); fmpz_print(N); flint_printf("\n");
            gr_mat_print(A, ctx); gr_mat_print(B, ctx);
            gr_mat_print(C, ctx); gr_mat_print(D, ctx);
            TEST_FUNCTION_FAIL("batched_mont != classical\n");
        }

        /* aliasing C == A, and squaring A*A */
        GR_MUST_SUCCEED(gr_mat_set(C, A, ctx));
        GR_MUST_SUCCEED(mpn_mod_mat_mul_batched_mont(C, C, B, ctx));
        if (gr_mat_equal(C, D, ctx) != T_TRUE)
            TEST_FUNCTION_FAIL("batched_mont aliasing (C==A) mismatch\n");

        GR_MUST_SUCCEED(gr_mat_mul_classical(D, A, A, ctx));
        GR_MUST_SUCCEED(mpn_mod_mat_mul_batched_mont(C, A, A, ctx));
        if (gr_mat_equal(C, D, ctx) != T_TRUE)
            TEST_FUNCTION_FAIL("batched_mont squaring mismatch\n");

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(C, ctx);
        gr_mat_clear(D, ctx);
        gr_ctx_clear(ctx);
        fmpz_clear(N);
    }

    TEST_FUNCTION_END(state);
}
