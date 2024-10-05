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

    for (iter = 0; iter < iters; iter++)
    {
        gr_poly_t A, B, Q, R, QBR;
        int status = GR_SUCCESS;
        gr_ctx_t ctx2;
        gr_ctx_struct * ctxptr;

        if (ctx == NULL)
        {
            gr_ctx_init_random(ctx2, state);
            ctxptr = ctx2;
        }
        else
            ctxptr = ctx;

        gr_poly_init(A, ctxptr);
        gr_poly_init(B, ctxptr);
        gr_poly_init(Q, ctxptr);
        gr_poly_init(R, ctxptr);
        gr_poly_init(QBR, ctxptr);

        status = GR_SUCCESS;

        status |= gr_poly_randtest(A, state, 1 + n_randint(state, maxn), ctxptr);
        status |= gr_poly_randtest(B, state, 1 + n_randint(state, maxn), ctxptr);
        if (A->length < B->length)
            gr_poly_swap(A, B, ctxptr);

        status |= gr_poly_randtest(Q, state, 1 + n_randint(state, maxn), ctxptr);
        status |= gr_poly_randtest(R, state, 1 + n_randint(state, maxn), ctxptr);

        /* randomly generate monic polynomials */
        if (n_randint(state, 2) && B->length >= 1)
            status |= gr_poly_set_coeff_si(B, B->length - 1, 1, ctxptr);

        if (n_randint(state, 3) == 0)
        {
            status |= gr_poly_mul(A, A, B, ctxptr);
            status |= gr_poly_add(A, A, R, ctxptr);
        }

        if (B->length >= 1)
        {
            gr_poly_fit_length(Q, A->length - B->length + 1, ctxptr);
            status |= div_impl(Q->coeffs, A->coeffs, A->length, B->coeffs, B->length, ctxptr);
            _gr_poly_set_length(Q, A->length - B->length + 1, ctxptr);
            _gr_poly_normalise(Q, ctxptr);

            status |= gr_poly_mul(R, Q, B, ctxptr);
            status |= gr_poly_sub(R, A, R, ctxptr);
        }
        else
        {
            status = GR_UNABLE;
        }

        if (status == GR_SUCCESS)
        {
            status |= gr_poly_mul(QBR, Q, B, ctxptr);
            status |= gr_poly_add(QBR, QBR, R, ctxptr);

            if (status == GR_SUCCESS && gr_poly_equal(QBR, A, ctxptr) == T_FALSE)
            {
                flint_printf("FAIL\n\n");
                flint_printf("A = "); gr_poly_print(A, ctxptr); flint_printf("\n");
                flint_printf("B = "); gr_poly_print(B, ctxptr); flint_printf("\n");
                flint_printf("Q = "); gr_poly_print(Q, ctxptr); flint_printf("\n");
                flint_printf("R = "); gr_poly_print(R, ctxptr); flint_printf("\n");
                flint_printf("Q*B + R = "); gr_poly_print(QBR, ctxptr); flint_printf("\n");
                flint_abort();
            }
        }

        gr_poly_clear(A, ctxptr);
        gr_poly_clear(B, ctxptr);
        gr_poly_clear(Q, ctxptr);
        gr_poly_clear(R, ctxptr);
        gr_poly_clear(QBR, ctxptr);

        if (ctx == NULL)
            gr_ctx_clear(ctxptr);
    }
}
