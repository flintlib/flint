/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_mat.h"
#include "gr_poly.h"
#include "fmpz_mat.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

static int
_my_polynomial_jet(gr_ptr res, gr_srcptr x, slong len, gr_poly_t f, gr_ctx_t ctx)
{
    slong i;
    int status = GR_SUCCESS;
    gr_poly_t g;
    ulong c = 1;

    gr_poly_init(g, ctx);
    status |= gr_poly_set(g, f, ctx);

    for (i = 0; i < len; i++)
    {
        c *= FLINT_MAX(1, i);
        status |= gr_poly_evaluate(GR_ENTRY(res, i, ctx->sizeof_elem), g, x, ctx);
        status |= gr_div_ui(GR_ENTRY(res, i, ctx->sizeof_elem),
                            GR_ENTRY(res, i, ctx->sizeof_elem), c, ctx);
        status |= gr_poly_derivative(g, g, ctx);
    }

    gr_poly_clear(g, ctx);
    return status;
}

TEST_GR_FUNCTION_START(gr_mat_func_jordan, state, count_success, count_domain, count_unable)
{
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong n;
        gr_ctx_t ctx;
        gr_mat_t A, B, C;
        gr_poly_t f;
        int status = GR_SUCCESS;

        while (1)
        {
            gr_ctx_init_random(ctx, state);
            if (gr_ctx_is_field(ctx) == T_TRUE)
                break;
            gr_ctx_clear(ctx);
        }

        if (ctx->methods == _ca_methods)
            n = n_randint(state, 3);
        else
            n = n_randint(state, 8);

        gr_mat_init(A, n, n, ctx);
        gr_mat_init(B, n, n, ctx);
        gr_mat_init(C, n, n, ctx);

        gr_poly_init(f, ctx);

        GR_MUST_SUCCEED(gr_poly_randtest(f, state, 1 + n_randint(state, 5), ctx));

        if (n_randint(state, 2))
        {
            GR_MUST_SUCCEED(gr_mat_randtest(A, state, ctx));
        }
        else
        {
            fmpq_mat_t Q;
            fmpq_mat_init(Q, n, n);
            fmpq_mat_randtest(Q, state, 5);
            status = gr_mat_set_fmpq_mat(A, Q, ctx);
            fmpq_mat_clear(Q);
        }

        status |= gr_mat_func_param_jordan(B, A, (gr_method_vec_scalar_op) _my_polynomial_jet, (gr_srcptr) f, ctx);

        if (status == GR_SUCCESS)
            status |= gr_mat_gr_poly_evaluate(C, f, A, ctx);

        if (status == GR_SUCCESS && gr_mat_equal(B, C, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n");
            gr_ctx_println(ctx);
            flint_printf("f = "), gr_poly_print(f, ctx); flint_printf("\n");
            flint_printf("A = "), gr_mat_print(A, ctx); flint_printf("\n");
            flint_printf("B = "), gr_mat_print(B, ctx); flint_printf("\n");
            flint_printf("C = "), gr_mat_print(C, ctx); flint_printf("\n");
            flint_abort();
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(C, ctx);

        gr_poly_clear(f, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
