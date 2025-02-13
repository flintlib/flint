/*
    Copyright (C) 2025 Lars Göttgens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "gr_mat.h"

TEST_FUNCTION_START(gr_mat_pow, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        slong n, i;
        ulong e_ui;
        fmpz_t e_fmpz, neg_e_fmpz;
        gr_ctx_t ctx;
        gr_mat_t A, B, C, Ainv;

        while (1)
        {
            gr_ctx_init_random(ctx, state);
            if (gr_ctx_is_finite(ctx) == T_TRUE)
                break;
            gr_ctx_clear(ctx);
        }

        n = n_randint(state, 5);
        e_ui = n_randint(state, 5);

        gr_mat_init(A, n, n, ctx);
        gr_mat_init(B, n, n, ctx);
        gr_mat_init(C, n, n, ctx);
        gr_mat_init(Ainv, n, n, ctx);

        GR_MUST_SUCCEED(gr_mat_randtest(A, state, ctx));
        GR_MUST_SUCCEED(gr_mat_randtest(B, state, ctx));

        GR_MUST_SUCCEED(gr_mat_pow_ui(B, A, e_ui, ctx));

        GR_MUST_SUCCEED(gr_mat_one(C, ctx));
        for (i = 0; i < e_ui; i++)
            GR_MUST_SUCCEED(gr_mat_mul(C, C, A, ctx));

        if (gr_mat_equal(C, B, ctx) == T_FALSE)
        {
            flint_printf("FAIL: results not equal (ui)\n");
            flint_printf("e %d", e_ui);
            gr_mat_print(A, ctx);
            gr_mat_print(C, ctx);
            gr_mat_print(B, ctx);
            fflush(stdout);
            flint_abort();
        }

        GR_MUST_SUCCEED(gr_mat_pow_ui(A, A, e_ui, ctx));

        if (gr_mat_equal(A, B, ctx) == T_FALSE)
        {
            flint_printf("FAIL: aliasing failed (ui)\n");
            fflush(stdout);
            flint_abort();
        }
        
        fmpz_init(e_fmpz);
        fmpz_init(neg_e_fmpz);
        fmpz_randtest(e_fmpz, state, 100);

        GR_MUST_SUCCEED(gr_mat_randtest(A, state, ctx));
        GR_MUST_SUCCEED(gr_mat_randtest(B, state, ctx));

        if (fmpz_is_zero(e_fmpz) || fmpz_sgn(e_fmpz) == 1)
        {
            GR_MUST_SUCCEED(gr_mat_pow_fmpz(B, A, e_fmpz, ctx));

            if (fmpz_fits_si(e_fmpz))
            {
                e_ui = fmpz_get_ui(e_fmpz);
                GR_MUST_SUCCEED(gr_mat_pow_ui(C, A, e_ui, ctx));

                if (gr_mat_equal(C, B, ctx) == T_FALSE)
                {
                    flint_printf("FAIL: results not equal (fmpz)\n");
                    fflush(stdout);
                    flint_abort();
                }
            }

            GR_MUST_SUCCEED(gr_mat_pow_fmpz(A, A, e_fmpz, ctx));

            if (gr_mat_equal(A, B, ctx) == T_FALSE)
            {
                flint_printf("FAIL: aliasing failed (fmpz)\n");
                fflush(stdout);
                flint_abort();
            }
        }
        else
        {
            int status_pow, status_inv;

            status_pow = gr_mat_pow_fmpz(B, A, e_fmpz, ctx);
            status_inv = gr_mat_inv(Ainv, A, ctx);
            if (status_pow != status_inv)
            {
                flint_printf("FAIL: pow with neg exp different from inv\n");
                fflush(stdout);
                flint_abort();
            }
            if (status_pow == GR_SUCCESS)
            {
                fmpz_neg(neg_e_fmpz, e_fmpz);
                GR_MUST_SUCCEED(gr_mat_pow_fmpz(C, Ainv, neg_e_fmpz, ctx));
                
                if (gr_mat_equal(C, B, ctx) == T_FALSE)
                {
                    flint_printf("FAIL: results not equal (neg fmpz)\n");
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        fmpz_clear(e_fmpz);
        fmpz_clear(neg_e_fmpz);
        
        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(C, ctx);
        gr_mat_clear(Ainv, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
