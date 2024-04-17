/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

void _gr_poly_test_div(gr_method_poly_binary_op div_impl,
    flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
{
    slong iter;
    gr_ctx_ptr given_ctx = ctx;

    for (iter = 0; iter < iters; iter++)
    {
        gr_poly_t A, B, Q, R, QBR;
        int status = GR_SUCCESS;
        gr_ctx_t my_ctx;
        gr_ctx_struct * ctx;

        if (given_ctx == NULL)
        {
            gr_ctx_init_random(my_ctx, state);
            ctx = my_ctx;
        }
        else
            ctx = given_ctx;

        gr_poly_init(A, ctx);
        gr_poly_init(B, ctx);
        gr_poly_init(Q, ctx);
        gr_poly_init(R, ctx);
        gr_poly_init(QBR, ctx);

        status = GR_SUCCESS;

        status |= gr_poly_randtest(A, state, 1 + n_randint(state, maxn), ctx);
        status |= gr_poly_randtest(B, state, 1 + n_randint(state, maxn), ctx);
        if (A->length < B->length)
            gr_poly_swap(A, B, ctx);

        status |= gr_poly_randtest(Q, state, 1 + n_randint(state, maxn), ctx);
        status |= gr_poly_randtest(R, state, 1 + n_randint(state, maxn), ctx);

        /* randomly generate monic polynomials */
        if (n_randint(state, 2) && B->length >= 1)
            status |= gr_poly_set_coeff_si(B, B->length - 1, 1, ctx);

        if (n_randint(state, 3) == 0)
        {
            status |= gr_poly_mul(A, A, B, ctx);
            status |= gr_poly_add(A, A, R, ctx);
        }

        if (B->length >= 1)
        {
            gr_poly_fit_length(Q, A->length - B->length + 1, ctx);
            status |= div_impl(Q->coeffs, A->coeffs, A->length, B->coeffs, B->length, ctx);
            _gr_poly_set_length(Q, A->length - B->length + 1, ctx);
            _gr_poly_normalise(Q, ctx);

            status |= gr_poly_mul(R, Q, B, ctx);
            status |= gr_poly_sub(R, A, R, ctx);
        }
        else
        {
            status = GR_UNABLE;
        }

        if (status == GR_SUCCESS)
        {
            status |= gr_poly_mul(QBR, Q, B, ctx);
            status |= gr_poly_add(QBR, QBR, R, ctx);

            if (status == GR_SUCCESS && gr_poly_equal(QBR, A, ctx) == T_FALSE)
            {
                flint_printf("FAIL\n\n");
                flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
                flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
                flint_printf("Q = "); gr_poly_print(Q, ctx); flint_printf("\n");
                flint_printf("R = "); gr_poly_print(R, ctx); flint_printf("\n");
                flint_printf("Q*B + R = "); gr_poly_print(QBR, ctx); flint_printf("\n");
                flint_abort();
            }
        }

        gr_poly_clear(A, ctx);
        gr_poly_clear(B, ctx);
        gr_poly_clear(Q, ctx);
        gr_poly_clear(R, ctx);
        gr_poly_clear(QBR, ctx);

        if (given_ctx == NULL)
            gr_ctx_clear(ctx);
    }
}
