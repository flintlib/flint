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
#include "util_test.h"

int
test_mul(flint_rand_t state, int which)
{
    gr_ctx_t ctx;
    slong n;
    gr_poly_t A, B, C, D;
    int status = GR_SUCCESS;

    gr_ctx_init_random(ctx, state);

    gr_poly_init(A, ctx);
    gr_poly_init(B, ctx);
    gr_poly_init(C, ctx);
    gr_poly_init(D, ctx);

    if (ctx->methods == _ca_methods)
        n = 4;
    else
        n = 10;

    GR_MUST_SUCCEED(gr_poly_randtest(A, state, 1 + n_randint(state, n), ctx));
    GR_MUST_SUCCEED(gr_poly_randtest(B, state, 1 + n_randint(state, n), ctx));
    GR_MUST_SUCCEED(gr_poly_randtest(C, state, 1 + n_randint(state, n), ctx));

    switch (which)
    {
        case 0:
            status |= gr_poly_mul_karatsuba(C, A, B, ctx);
            break;
        case 1:
            status |= gr_poly_set(C, A, ctx);
            status |= gr_poly_mul_karatsuba(C, C, B, ctx);
            break;
        case 2:
            status |= gr_poly_set(C, B, ctx);
            status |= gr_poly_mul_karatsuba(C, A, C, ctx);
            break;
        case 3:
            status |= gr_poly_set(B, A, ctx);
            status |= gr_poly_mul_karatsuba(C, A, A, ctx);
            break;
        case 4:
            status |= gr_poly_set(B, A, ctx);
            status |= gr_poly_set(C, A, ctx);
            status |= gr_poly_mul_karatsuba(C, C, C, ctx);
            break;

        default:
            flint_abort();
    }

    /* todo: should explicitly call basecase mul */
    status |= gr_poly_mullow(D, A, B, FLINT_MAX(0, A->length + B->length - 1), ctx);

    if (status == GR_SUCCESS && gr_poly_equal(C, D, ctx) == T_FALSE)
    {
        flint_printf("FAIL\n\n");
        flint_printf("which = %d, n = %wd\n\n", which, n);
        gr_ctx_println(ctx);
        flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n\n");
        flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n\n");
        flint_printf("C = "); gr_poly_print(C, ctx); flint_printf("\n\n");
        flint_printf("D = "); gr_poly_print(D, ctx); flint_printf("\n\n");
        flint_abort();
    }

    gr_poly_clear(A, ctx);
    gr_poly_clear(B, ctx);
    gr_poly_clear(C, ctx);
    gr_poly_clear(D, ctx);

    gr_ctx_clear(ctx);

    return status;
}

TEST_FUNCTION_START(gr_poly_mul_karatsuba, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        test_mul(state, n_randint(state, 5));
    }

    TEST_FUNCTION_END(state);
}
