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

TEST_FUNCTION_START(mpn_mod_mat_pow_ui_batched_mont, state)
{
    slong iter;

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        fmpz_t N, val;
        gr_mat_t A, C, D;
        slong i, j, bits;
        ulong exp;

        /* odd modulus with 65 <= bits <= 248 (>= 2 limbs, < 2^248) */
        bits = 65 + n_randint(state, 248 - 65 + 1);
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
        gr_mat_init(C, 4, 4, ctx);
        gr_mat_init(D, 4, 4, ctx);

        fmpz_init(val);
        for (i = 0; i < 4; i++)
            for (j = 0; j < 4; j++)
            {
                fmpz_randm(val, state, N);
                GR_MUST_SUCCEED(gr_set_fmpz(gr_mat_entry_ptr(A, i, j, ctx), val, ctx));
            }
        fmpz_clear(val);

        /* mix of edge and random exponents, including 0 and 1 */
        switch (n_randint(state, 4))
        {
            case 0:  exp = 0; break;
            case 1:  exp = 1; break;
            case 2:  exp = n_randint(state, 64); break;
            default: exp = n_randtest(state); break;   /* full-width chain */
        }

        GR_MUST_SUCCEED(mpn_mod_mat_pow_ui_batched_mont(C, A, exp, ctx));
        GR_MUST_SUCCEED(gr_mat_pow_ui(D, A, exp, ctx));

        if (gr_mat_equal(C, D, ctx) != T_TRUE)
        {
            flint_printf("FAIL: iter %wd, exp = %wu\n", iter, exp);
            flint_printf("n = "); fmpz_print(N); flint_printf("\n");
            gr_mat_print(A, ctx);
            gr_mat_print(C, ctx);
            gr_mat_print(D, ctx);
            TEST_FUNCTION_FAIL("batched-mont A^exp != gr_mat_pow_ui\n");
        }

        /* aliasing: res == A */
        GR_MUST_SUCCEED(gr_mat_set(C, A, ctx));
        GR_MUST_SUCCEED(mpn_mod_mat_pow_ui_batched_mont(C, C, exp, ctx));
        if (gr_mat_equal(C, D, ctx) != T_TRUE)
        {
            flint_printf("FAIL (aliased res==A): iter %wd, exp = %wu\n", iter, exp);
            flint_printf("n = "); fmpz_print(N); flint_printf("\n");
            TEST_FUNCTION_FAIL("batched-mont aliasing mismatch\n");
        }

        gr_mat_clear(A, ctx);
        gr_mat_clear(C, ctx);
        gr_mat_clear(D, ctx);
        gr_ctx_clear(ctx);
        fmpz_clear(N);
    }

    TEST_FUNCTION_END(state);
}
