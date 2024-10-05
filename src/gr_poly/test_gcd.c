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

static int
gr_poly_gcd_wrapper(gr_method_poly_gcd_op gcd_impl, gr_poly_t G, const gr_poly_t A,
                        const gr_poly_t B, gr_ctx_t ctx)
{
    slong lenA = A->length, lenB = B->length, lenG;
    slong sz = ctx->sizeof_elem;
    gr_ptr g;
    int status;

    if (A->length == 0 && B->length == 0)
        return gr_poly_zero(G, ctx);

    if (A->length == 0)
        return gr_poly_make_monic(G, B, ctx);

    if (B->length == 0)
        return gr_poly_make_monic(G, A, ctx);

    if (A->length < B->length)
        return gr_poly_gcd_wrapper(gcd_impl, G, B, A, ctx);

    if (gr_is_zero(GR_ENTRY(A->coeffs, A->length - 1, sz), ctx) != T_FALSE ||
        gr_is_zero(GR_ENTRY(B->coeffs, B->length - 1, sz), ctx) != T_FALSE)
    {
        return GR_UNABLE;
    }

    /* lenA >= lenB >= 1 */
    if (G == A || G == B)
    {
        g = flint_malloc(FLINT_MIN(lenA, lenB) * sz);
        _gr_vec_init(g, FLINT_MIN(lenA, lenB), ctx);
    }
    else
    {
        gr_poly_fit_length(G, FLINT_MIN(lenA, lenB), ctx);
        g = G->coeffs;
    }

    status = gcd_impl(g, &lenG, A->coeffs, lenA, B->coeffs, lenB, ctx);

    if (G == A || G == B)
    {
        _gr_vec_clear(G->coeffs, G->alloc, ctx);
        flint_free(G->coeffs);
        G->coeffs = g;
        G->alloc = FLINT_MIN(lenA, lenB);
        G->length = FLINT_MIN(lenA, lenB);
    }
    _gr_poly_set_length(G, lenG, ctx);

    if (status == GR_SUCCESS && lenG != 0)
    {
        if (lenG == 1)
            status = gr_one(G->coeffs, ctx);
        else
            status = gr_poly_make_monic(G, G, ctx);
    }

    return status;
}

void _gr_poly_test_gcd(gr_method_poly_gcd_op gcd_impl,
    flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
{
    slong iter;

    for (iter = 0; iter < iters; iter++)
    {
        gr_ctx_t ctx2;
        gr_ctx_struct * ctxptr;

        if (ctx == NULL)
        {
            gr_ctx_init_random(ctx2, state);
            ctxptr = ctx2;
        }
        else
            ctxptr = ctx;

        {
            gr_poly_t A, B, C, AC, BC, MC, G;
            int status = GR_SUCCESS;
            slong n;

            gr_poly_init(A, ctxptr);
            gr_poly_init(B, ctxptr);
            gr_poly_init(C, ctxptr);
            gr_poly_init(AC, ctxptr);
            gr_poly_init(BC, ctxptr);
            gr_poly_init(MC, ctxptr);
            gr_poly_init(G, ctxptr);

            n = 1 + n_randint(state, maxn);

            status = gr_poly_randtest(A, state, n, ctxptr);
            status |= gr_poly_randtest(B, state, n, ctxptr);

            status |= gr_poly_gcd_wrapper(gcd_impl, G, A, B, ctxptr);

            if (status == GR_SUCCESS)
            {
                if (gr_poly_is_zero(G, ctxptr) == T_FALSE)
                {
                    status |= gr_poly_divrem(AC, BC, A, G, ctxptr);

                    if (status == GR_SUCCESS && gr_poly_is_zero(BC, ctxptr) == T_FALSE)
                    {
                        flint_printf("FAIL: gcd does not divide A\n\n");
                        flint_printf("A = "); gr_poly_print(A, ctxptr); flint_printf("\n");
                        flint_printf("B = "); gr_poly_print(B, ctxptr); flint_printf("\n");
                        flint_printf("G = "); gr_poly_print(G, ctxptr); flint_printf("\n");
                        flint_abort();
                    }

                    status |= gr_poly_divrem(AC, BC, B, G, ctxptr);

                    if (status == GR_SUCCESS && gr_poly_is_zero(BC, ctxptr) == T_FALSE)
                    {
                        flint_printf("FAIL: gcd does not divide B\n\n");
                        flint_printf("A = "); gr_poly_print(A, ctxptr); flint_printf("\n");
                        flint_printf("B = "); gr_poly_print(B, ctxptr); flint_printf("\n");
                        flint_printf("G = "); gr_poly_print(G, ctxptr); flint_printf("\n");
                        flint_abort();
                    }
                }
            }

            if (status == GR_SUCCESS && gr_poly_is_one(G, ctxptr) == T_TRUE)
            {
                status |= gr_poly_randtest(C, state, n, ctxptr);

                status |= gr_poly_mul(AC, A, C, ctxptr);
                status |= gr_poly_mul(BC, B, C, ctxptr);

                switch (n_randint(state, 3))
                {
                    case 0:
                        status |= gr_poly_set(G, AC, ctxptr);
                        status |= gr_poly_gcd(G, G, BC, ctxptr);
                        break;
                    case 1:
                        status |= gr_poly_set(G, BC, ctxptr);
                        status |= gr_poly_gcd(G, AC, G, ctxptr);
                        break;
                    default:
                        status |= gr_poly_gcd(G, AC, BC, ctxptr);
                        break;
                }

                if (status == GR_SUCCESS)
                {
                    if (gr_poly_is_zero(C, ctxptr) == T_FALSE)
                    {
                        status |= gr_poly_make_monic(MC, C, ctxptr);

                        if (status == GR_SUCCESS && gr_poly_equal(MC, G, ctxptr) == T_FALSE)
                        {
                            flint_printf("FAIL\n\n");
                            flint_printf("A = "); gr_poly_print(A, ctxptr); flint_printf("\n");
                            flint_printf("B = "); gr_poly_print(B, ctxptr); flint_printf("\n");
                            flint_printf("C = "); gr_poly_print(C, ctxptr); flint_printf("\n");
                            flint_printf("G = "); gr_poly_print(G, ctxptr); flint_printf("\n");
                            flint_abort();
                        }
                    }
                }
            }

            status = gr_poly_randtest(A, state, n, ctxptr);
            status |= gr_poly_gcd_wrapper(gcd_impl, G, A, A, ctxptr);

            if (status == GR_SUCCESS)
            {
                status |= gr_poly_make_monic(B, A, ctxptr);

                if (status == GR_SUCCESS && gr_poly_equal(G, B, ctxptr) == T_FALSE)
                {
                    flint_printf("FAIL (self)\n\n");
                    flint_printf("A = "); gr_poly_print(A, ctxptr); flint_printf("\n");
                    flint_printf("B = "); gr_poly_print(B, ctxptr); flint_printf("\n");
                    flint_printf("G = "); gr_poly_print(G, ctxptr); flint_printf("\n");
                    flint_abort();
                }
            }

            gr_poly_clear(A, ctxptr);
            gr_poly_clear(B, ctxptr);
            gr_poly_clear(C, ctxptr);
            gr_poly_clear(AC, ctxptr);
            gr_poly_clear(BC, ctxptr);
            gr_poly_clear(MC, ctxptr);
            gr_poly_clear(G, ctxptr);
        }

        if (ctx == NULL)
            gr_ctx_clear(ctxptr);
    }
}
