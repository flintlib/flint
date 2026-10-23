/*
    Copyright (C) 2026 Maria Neagoie

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_ore_poly.h"

TEST_GR_FUNCTION_START(gr_ore_poly_divrem, state, count_success, count_domain, count_unable)
{
    slong iter;

    /* regression test: the remainder must be fully reduced */
    {
        int status = GR_SUCCESS;
        gr_ctx_t QQ;
        gr_ore_poly_ctx_t ore_ctx;
        gr_ore_poly_t U, V, Q, R, Q2, R2;

        gr_ctx_init_fmpq(QQ);
        gr_ore_poly_ctx_init(ore_ctx, QQ, 0, ORE_ALGEBRA_COMMUTATIVE);
        gr_ore_poly_init(U, ore_ctx);
        gr_ore_poly_init(V, ore_ctx);
        gr_ore_poly_init(Q, ore_ctx);
        gr_ore_poly_init(R, ore_ctx);
        gr_ore_poly_init(Q2, ore_ctx);
        gr_ore_poly_init(R2, ore_ctx);

        status |= gr_ore_poly_set_str(U, "D^2 + D + 1", ore_ctx);
        status |= gr_ore_poly_set_str(V, "D", ore_ctx);
        status |= gr_ore_poly_set_str(Q2, "D + 1", ore_ctx);
        status |= gr_ore_poly_one(R2, ore_ctx);
        status |= gr_ore_poly_divrem(Q, R, U, V, ore_ctx);

        if (status != GR_SUCCESS || gr_ore_poly_equal(Q, Q2, ore_ctx) != T_TRUE ||
                gr_ore_poly_equal(R, R2, ore_ctx) != T_TRUE)
        {
            flint_printf("FAIL (D^2 + D + 1) / D\n\n");
            flint_printf("Q = "); gr_ore_poly_print(Q, ore_ctx); flint_printf("\n");
            flint_printf("R = "); gr_ore_poly_print(R, ore_ctx); flint_printf("\n");
            flint_abort();
        }

        gr_ore_poly_clear(U, ore_ctx);
        gr_ore_poly_clear(V, ore_ctx);
        gr_ore_poly_clear(Q, ore_ctx);
        gr_ore_poly_clear(R, ore_ctx);
        gr_ore_poly_clear(Q2, ore_ctx);
        gr_ore_poly_clear(R2, ore_ctx);
        gr_ore_poly_ctx_clear(ore_ctx);
        gr_ctx_clear(QQ);
    }

    for (iter = 0; iter < 3000 * flint_test_multiplier(); iter++)
    {
        int status = GR_SUCCESS;
        gr_ctx_t ctx;
        gr_ore_poly_ctx_t ore_ctx;
        gr_ore_poly_t A, B, Q, R, QBR;

        gr_ore_poly_ctx_init_randtest2(ctx, ore_ctx, state);

        gr_ore_poly_init(A, ore_ctx);
        gr_ore_poly_init(B, ore_ctx);
        gr_ore_poly_init(Q, ore_ctx);
        gr_ore_poly_init(R, ore_ctx);
        gr_ore_poly_init(QBR, ore_ctx);

        status |= gr_ore_poly_randtest(A, state, 1 + n_randint(state, 6), ore_ctx);
        status |= gr_ore_poly_randtest(B, state, 1 + n_randint(state, 6), ore_ctx);
        status |= gr_ore_poly_randtest(Q, state, 1 + n_randint(state, 6), ore_ctx);
        status |= gr_ore_poly_randtest(R, state, 1 + n_randint(state, 6), ore_ctx);

        if (n_randint(state, 3) == 0)
        {
            status |= gr_ore_poly_mul(A, A, B, ore_ctx);
            status |= gr_ore_poly_add(A, A, R, ore_ctx);
        }

        if (status == GR_SUCCESS)
        {
            // test aliasing
            switch (n_randint(state, 5))
            {
                case 0:
                    status |= gr_ore_poly_set(Q, A, ore_ctx);
                    status |= gr_ore_poly_divrem(Q, R, Q, B, ore_ctx);
                    break;
                case 1:
                    status |= gr_ore_poly_set(R, A, ore_ctx);
                    status |= gr_ore_poly_divrem(Q, R, R, B, ore_ctx);
                    break;
                case 2:
                    status |= gr_ore_poly_set(Q, B, ore_ctx);
                    status |= gr_ore_poly_divrem(Q, R, A, Q, ore_ctx);
                    break;
                case 3:
                    status |= gr_ore_poly_set(R, B, ore_ctx);
                    status |= gr_ore_poly_divrem(Q, R, A, R, ore_ctx);
                    break;
                default:
                    status |= gr_ore_poly_divrem(Q, R, A, B, ore_ctx);
                    break;
            }

            if (status == GR_SUCCESS)
            {
                status |= gr_ore_poly_mul(QBR, Q, B, ore_ctx);
                status |= gr_ore_poly_add(QBR, QBR, R, ore_ctx);

                if (status == GR_SUCCESS && (gr_ore_poly_equal(QBR, A, ore_ctx) == T_FALSE ||
                        gr_ore_poly_length(R, ore_ctx) >= FLINT_MAX(gr_ore_poly_length(B, ore_ctx), 1)))
                {
                    flint_printf("FAIL\n\n");
                    flint_printf("A = "); gr_ore_poly_print(A, ore_ctx); flint_printf("\n");
                    flint_printf("B = "); gr_ore_poly_print(B, ore_ctx); flint_printf("\n");
                    flint_printf("Q = "); gr_ore_poly_print(Q, ore_ctx); flint_printf("\n");
                    flint_printf("R = "); gr_ore_poly_print(R, ore_ctx); flint_printf("\n");
                    flint_printf("Q*B + R = "); gr_ore_poly_print(QBR, ore_ctx); flint_printf("\n");
                    flint_abort();
                }
            }
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_ore_poly_clear(A, ore_ctx);
        gr_ore_poly_clear(B, ore_ctx);
        gr_ore_poly_clear(Q, ore_ctx);
        gr_ore_poly_clear(R, ore_ctx);
        gr_ore_poly_clear(QBR, ore_ctx);

        gr_ore_poly_ctx_clear(ore_ctx);
        gr_ctx_clear(ctx);
    }
    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
