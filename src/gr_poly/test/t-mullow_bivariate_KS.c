/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_poly.h"

TEST_FUNCTION_START(gr_poly_mullow_bivariate_KS, state)
{
    slong iter;

    for (iter = 0; iter < 10000; iter++)
    {
        gr_ctx_t ctx, ctx2;
        slong n;
        gr_poly_t A, B, C, D;
        int status = GR_SUCCESS;
        int aliasing = n_randint(state, 5);
        int karatsuba = n_randint(state, 2);

        while (1)
        {
            gr_ctx_init_random(ctx2, state);
            if (ctx2->methods != _ca_methods)
                break;
            gr_ctx_clear(ctx2);
        }

        if (n_randint(state, 2))
        {
            gr_ctx_init_gr_poly(ctx, ctx2);
        }
        else
        {
            FLINT_SWAP(gr_ctx_struct, *ctx, *ctx2);
            gr_ctx_init_fmpz(ctx2); /* dummy */
        }

        gr_poly_init(A, ctx);
        gr_poly_init(B, ctx);
        gr_poly_init(C, ctx);
        gr_poly_init(D, ctx);

        n = n_randint(state, 6);

        GR_MUST_SUCCEED(gr_poly_randtest(A, state, 1 + n_randint(state, 6), ctx));
        GR_MUST_SUCCEED(gr_poly_randtest(B, state, 1 + n_randint(state, 6), ctx));
        GR_MUST_SUCCEED(gr_poly_randtest(C, state, 1 + n_randint(state, 6), ctx));
        GR_MUST_SUCCEED(gr_poly_randtest(D, state, 1 + n_randint(state, 6), ctx));

        switch (aliasing)
        {
            case 0:
                status |= gr_poly_mullow_bivariate_KS(C, A, B, n, ctx);
                break;
            case 1:
                status |= gr_poly_set(C, A, ctx);
                status |= gr_poly_mullow_bivariate_KS(C, C, B, n, ctx);
                break;
            case 2:
                status |= gr_poly_set(C, B, ctx);
                status |= gr_poly_mullow_bivariate_KS(C, A, C, n, ctx);
                break;
            case 3:
                status |= gr_poly_set(B, A, ctx);
                status |= gr_poly_mullow_bivariate_KS(C, A, A, n, ctx);
                break;
            case 4:
                status |= gr_poly_set(B, A, ctx);
                status |= gr_poly_set(C, A, ctx);
                status |= gr_poly_mullow_bivariate_KS(C, C, C, n, ctx);
                break;

            default:
                flint_abort();
        }

        if (status == GR_SUCCESS)
        {
            status |= gr_poly_mullow_classical(D, A, B, n, ctx);

            if (status == GR_SUCCESS && gr_poly_equal(C, D, ctx) == T_FALSE)
            {
                flint_printf("FAIL\n\n");
                flint_printf("aliasing = %d, n = %wd, karatsuba = %d\n\n", aliasing, n, karatsuba);
                gr_ctx_println(ctx);
                flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n\n");
                flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n\n");
                flint_printf("C = "); gr_poly_print(C, ctx); flint_printf("\n\n");
                flint_printf("D = "); gr_poly_print(D, ctx); flint_printf("\n\n");
                flint_abort();
            }
        }

        gr_poly_clear(A, ctx);
        gr_poly_clear(B, ctx);
        gr_poly_clear(C, ctx);
        gr_poly_clear(D, ctx);

        gr_ctx_clear(ctx);
        gr_ctx_clear(ctx2);
    }

    TEST_FUNCTION_END(state);
}

