/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_poly.h"

TEST_FUNCTION_START(gr_poly_divexact, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        int status;
        gr_ctx_t ctx;
        gr_poly_t A, B, Q, Q2;

        gr_ctx_init_random(ctx, state);

        gr_poly_init(A, ctx);
        gr_poly_init(B, ctx);
        gr_poly_init(Q, ctx);
        gr_poly_init(Q2, ctx);

        status = GR_SUCCESS;

        status |= gr_poly_randtest(B, state, 1 + n_randint(state, 6), ctx);

        if (gr_poly_is_zero(B, ctx) == T_FALSE)
        {
            status |= gr_poly_randtest(Q, state, 1 + n_randint(state, 6), ctx);
            status |= gr_poly_randtest(Q2, state, 1 + n_randint(state, 6), ctx);
            status |= gr_poly_mul(A, Q2, B, ctx);

            switch (n_randint(state, 12))
            {
                case 0:
                    status |= gr_poly_set(Q, A, ctx);
                    status |= gr_poly_divexact(Q, Q, B, ctx);
                    break;
                case 1:
                    status |= gr_poly_set(Q, B, ctx);
                    status |= gr_poly_divexact(Q, A, B, ctx);
                    break;
                case 2:
                    status |= gr_poly_divexact(Q, A, B, ctx);
                    break;

                case 3:
                    status |= gr_poly_set(Q, A, ctx);
                    status |= gr_poly_divexact_basecase(Q, Q, B, ctx);
                    break;
                case 4:
                    status |= gr_poly_set(Q, B, ctx);
                    status |= gr_poly_divexact_basecase(Q, A, B, ctx);
                    break;
                case 5:
                    status |= gr_poly_divexact_basecase(Q, A, B, ctx);
                    break;

                case 6:
                    status |= gr_poly_set(Q, A, ctx);
                    status |= gr_poly_divexact_bidirectional(Q, Q, B, ctx);
                    break;
                case 7:
                    status |= gr_poly_set(Q, B, ctx);
                    status |= gr_poly_divexact_bidirectional(Q, A, B, ctx);
                    break;
                case 8:
                    status |= gr_poly_divexact_bidirectional(Q, A, B, ctx);
                    break;

                case 9:
                    status |= gr_poly_set(Q, A, ctx);
                    status |= gr_poly_divexact_basecase_bidirectional(Q, Q, B, ctx);
                    break;
                case 10:
                    status |= gr_poly_set(Q, B, ctx);
                    status |= gr_poly_divexact_basecase_bidirectional(Q, A, B, ctx);
                    break;
                default:
                    status |= gr_poly_divexact_basecase_bidirectional(Q, A, B, ctx);
                    break;
            }

            if (gr_ctx_is_integral_domain(ctx) == T_TRUE &&
                ((status == GR_SUCCESS && gr_poly_equal(Q, Q2, ctx) == T_FALSE)
                    || status == GR_DOMAIN
                    || (ctx->which_ring == GR_CTX_FMPZ && status != GR_SUCCESS)))
            {
                flint_printf("FAIL\n\n");
                flint_printf("%d\n", status);
                gr_ctx_println(ctx);
                flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
                flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
                flint_printf("Q = "); gr_poly_print(Q, ctx); flint_printf("\n");
                flint_printf("Q2 = "); gr_poly_print(Q2, ctx); flint_printf("\n");
                flint_abort();
            }
        }

        gr_poly_clear(A, ctx);
        gr_poly_clear(B, ctx);
        gr_poly_clear(Q, ctx);
        gr_poly_clear(Q2, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
