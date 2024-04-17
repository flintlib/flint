/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/


#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_xgcd_generic(slong * lenG, gr_ptr G, gr_ptr S, gr_ptr T, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
{
    FLINT_ASSERT(lenA >= lenB);
    FLINT_ASSERT(lenB >= 1);

    return _gr_poly_xgcd_euclidean(lenG, G, S, T, A, lenA, B, lenB, ctx);
}

int
gr_poly_xgcd_wrapper(gr_method_poly_xgcd_op _xgcd_op, gr_poly_t G, gr_poly_t S, gr_poly_t T, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
{
    if (A->length < B->length)
    {
        return gr_poly_xgcd_wrapper(_xgcd_op, G, T, S, B, A, ctx);
    }
    else                        /* lenA >= lenB >= 0 */
    {
        const slong lenA = A->length, lenB = B->length;
        int status = GR_SUCCESS;
        slong sz = ctx->sizeof_elem;
        gr_ptr t;


        if (lenA == 0)          /* lenA = lenB = 0 */
        {
            status |= gr_poly_zero(G, ctx);
            status |= gr_poly_zero(S, ctx);
            status |= gr_poly_zero(T, ctx);
        }
        else if (lenB == 0)     /* lenA > lenB = 0 */
        {
            GR_TMP_INIT(t, ctx);
            status |= gr_inv(t, GR_ENTRY(A->coeffs, lenA - 1, sz), ctx);
            status |= gr_poly_mul_scalar(G, A, t, ctx);
            status |= gr_poly_zero(T, ctx);
            status |= gr_poly_set_scalar(S, t, ctx);
            GR_TMP_CLEAR(t, ctx);
        }
        else if (gr_is_zero(GR_ENTRY(A->coeffs, A->length - 1, sz), ctx) != T_FALSE ||
                gr_is_zero(GR_ENTRY(B->coeffs, B->length - 1, sz), ctx) != T_FALSE)
        {
            status |= GR_UNABLE;
        }
        else if (lenB == 1)  /* lenA >= lenB = 1 */
        {
            GR_TMP_INIT(t, ctx);
            status |= gr_inv(t, B->coeffs, ctx);
            status |= gr_poly_set_scalar(T, t, ctx);
            status |= gr_poly_one(G, ctx);
            status |= gr_poly_zero(S, ctx);
            GR_TMP_CLEAR(t, ctx);
        }
        else                    /* lenA >= lenB >= 2 */
        {
            gr_ptr g, *s, *t;
            slong lenG;

            if (G == A || G == B)
            {
                g = flint_malloc(sz * FLINT_MIN(lenA, lenB));
                _gr_vec_init(g, FLINT_MIN(lenA, lenB), ctx);
            }
            else
            {
                gr_poly_fit_length(G, FLINT_MIN(lenA, lenB), ctx);
                g = G->coeffs;
            }

            if (S == A || S == B)
            {
                s = flint_malloc(sz * lenB);
                _gr_vec_init(s, lenB, ctx);
            }
            else
            {
                gr_poly_fit_length(S, lenB, ctx);
                s = S->coeffs;
            }

            if (T == A || T == B)
            {
                t = flint_malloc(sz * lenA);
                _gr_vec_init(t, lenA, ctx);
            }
            else
            {
                gr_poly_fit_length(T, lenA, ctx);
                t = T->coeffs;
            }

            status |= _xgcd_op(&lenG, g, s, t, A->coeffs, lenA, B->coeffs, lenB, ctx);

            if (G == A || G == B)
            {
                _gr_vec_clear(G->coeffs, G->alloc, ctx);
                flint_free(G->coeffs);
                G->coeffs = g;
                G->alloc = FLINT_MIN(lenA, lenB);
                G->length = G->alloc;
            }

            if (S == A || S == B)
            {
                _gr_vec_clear(S->coeffs, S->alloc, ctx);
                flint_free(S->coeffs);
                S->coeffs = s;
                S->alloc = lenB;
                S->length = S->alloc;
            }

            if (T == A || T == B)
            {
                _gr_vec_clear(T->coeffs, T->alloc, ctx);
                flint_free(T->coeffs);
                T->coeffs = t;
                T->alloc = lenA;
                T->length = T->alloc;
            }

            _gr_poly_set_length(G, lenG, ctx);
            _gr_poly_set_length(S, FLINT_MAX(lenB - lenG, 1), ctx);
            _gr_poly_set_length(T, FLINT_MAX(lenA - lenG, 1), ctx);
            _gr_poly_normalise(S, ctx);
            _gr_poly_normalise(T, ctx);

            if (status == GR_SUCCESS && gr_is_one(GR_ENTRY(G->coeffs, G->length - 1, sz), ctx) != T_TRUE)
            {
                GR_TMP_INIT(t, ctx);
                status |= gr_inv(t, GR_ENTRY(G->coeffs, G->length - 1, sz), ctx);
                status |= gr_poly_mul_scalar(G, G, t, ctx);
                status |= gr_poly_mul_scalar(S, S, t, ctx);
                status |= gr_poly_mul_scalar(T, T, t, ctx);
                GR_TMP_CLEAR(t, ctx);
            }
        }

        return status;
    }
}

int
gr_poly_xgcd(gr_poly_t G, gr_poly_t S, gr_poly_t T, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
{
    return gr_poly_xgcd_wrapper((gr_method_poly_xgcd_op) _gr_poly_xgcd, G, S, T, A, B, ctx);
}
