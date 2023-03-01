/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/


#include "gr_vec.h"
#include "gr_poly.h"

#define GR_VEC_NORM(status, R, lenR, sz, ctx) \
    do { \
        while ((lenR) > 0) \
        { \
            truth_t is_zero; \
            is_zero = gr_is_zero(GR_ENTRY((R), (lenR) - 1, sz), (ctx)); \
            if (is_zero == T_TRUE) \
                (lenR)--; \
            else if (is_zero == T_UNKNOWN) \
            { \
                (status) |= GR_UNABLE; \
                break; \
            } \
            else \
            { \
                break; \
            } \
        } \
    } while (0) \

#define GR_VEC_SWAP(vec1, len1, vec2, len2) \
do {                                          \
    gr_ptr __t;                                \
    slong __tn;                                \
    __t    = (vec1);                          \
    (vec1) = (vec2);                          \
    (vec2) = __t;                             \
    __tn   = (len1);                          \
    (len1) = (len2);                          \
    (len2) = __tn;                            \
} while (0);

/* todo: use invB */
int
_gr_poly_xgcd_euclidean(slong * lenG, gr_ptr G, gr_ptr S, gr_ptr T, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_srcptr invB, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    status |= _gr_vec_zero(G, lenB, ctx);
    status |= _gr_vec_zero(S, lenB - 1, ctx);
    status |= _gr_vec_zero(T, lenA - 1, ctx);

    if (lenB == 1)
    {
        status |= gr_set(G, B, ctx);
        status |= gr_one(T, ctx);
        *lenG = 1;
    }
    else
    {
        gr_ptr Q, R;
        slong lenQ, lenR;

        GR_TMP_INIT_VEC(Q, 2 * lenA, ctx);
        R = GR_ENTRY(Q, lenA, sz);

        /* todo: preinv1 B */
        status |= _gr_poly_divrem(Q, R, A, lenA, B, lenB, ctx);
        lenR = lenB - 1;
        GR_VEC_NORM(status, R, lenR, sz, ctx);

        if (lenR == 0)
        {
            status |= _gr_vec_set(G, B, lenB, ctx);
            status |= gr_one(T, ctx);
            GR_TMP_CLEAR_VEC(Q, 2 * lenA, ctx);
            *lenG = lenB;
            return status;
        }
        else
        {
            gr_ptr D, U, V1, V3, W;
            slong lenD, lenU, lenV1, lenV3, lenW;

            GR_TMP_INIT_VEC(W, FLINT_MAX(5 * lenB, lenA + lenB), ctx);
            D = GR_ENTRY(W, lenB, sz);
            U = GR_ENTRY(D, lenB, sz);
            V1 = GR_ENTRY(U, lenB, sz);
            V3 = GR_ENTRY(V1, lenB, sz);

            lenU = 0;
            status |= _gr_vec_set(D, B, lenB, ctx);
            lenD = lenB;
            status |= gr_one(V1, ctx);
            lenV1 = 1;
            lenV3 = 0;
            GR_VEC_SWAP(V3, lenV3, R, lenR);

            do
            {
                /* todo: in-place remainder (see fmpz_mod_poly) */
                status |= _gr_poly_divrem(Q, R, D, lenD, V3, lenV3, ctx);
                lenQ = lenD - lenV3 + 1;
                lenR = lenV3 - 1;
                GR_VEC_NORM(status, R, lenR, sz, ctx);                if (status != GR_SUCCESS)
                    break;

                if (lenV1 >= lenQ)
                    status |= _gr_poly_mul(W, V1, lenV1, Q, lenQ, ctx);
                else
                    status |= _gr_poly_mul(W, Q, lenQ, V1, lenV1, ctx);
                lenW = lenQ + lenV1 - 1;

                status |= _gr_poly_sub(U, U, lenU, W, lenW, ctx);
                lenU = FLINT_MAX(lenU, lenW);
                GR_VEC_NORM(status, U, lenU, sz, ctx);
                if (status != GR_SUCCESS)
                    break;

                GR_VEC_SWAP(U, lenU, V1, lenV1);
                {
                    gr_ptr __t;
                    slong __tn;

                    __t = D;
                    D = V3;
                    V3 = R;
                    R = __t;
                    __tn = lenD;
                    lenD = lenV3;
                    lenV3 = lenR;
                    lenR = __tn;
                }

            } while (lenV3 != 0);

            status |= _gr_vec_set(G, D, lenD, ctx);
            status |= _gr_vec_set(S, U, lenU, ctx);

            {
                lenQ = lenA + lenU - 1;

                status |= _gr_poly_mul(Q, A, lenA, S, lenU, ctx);
                status |= _gr_vec_neg(Q, Q, lenQ, ctx);
                status |= _gr_poly_add(Q, G, lenD, Q, lenQ, ctx);
                /* todo: preinv1 B */
                status |= _gr_poly_divrem(T, W, Q, lenQ, B, lenB, ctx);
            }

            GR_TMP_CLEAR_VEC(W, FLINT_MAX(5 * lenB, lenA + lenB), ctx);
            GR_TMP_CLEAR_VEC(Q, 2 * lenA, ctx);
            *lenG = lenD;
        }
    }

    return status;
}

int
gr_poly_xgcd_euclidean(gr_poly_t G, gr_poly_t S, gr_poly_t T, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
{
    if (A->length < B->length)
    {
        return gr_poly_xgcd_euclidean(G, T, S, B, A, ctx);
    }
    else                        /* lenA >= lenB >= 0 */
    {
        const slong lenA = A->length, lenB = B->length;
        int status = GR_SUCCESS;
        slong sz = ctx->sizeof_elem;
        gr_ptr inv;

        GR_TMP_INIT(inv, ctx);

        if (lenA == 0)          /* lenA = lenB = 0 */
        {
            status |= gr_poly_zero(G, ctx);
            status |= gr_poly_zero(S, ctx);
            status |= gr_poly_zero(T, ctx);
        }
        else if (lenB == 0)     /* lenA > lenB = 0 */
        {
            status |= gr_inv(inv, GR_ENTRY(A->coeffs, lenA - 1, sz), ctx);
            status |= gr_poly_mul_scalar(G, A, inv, ctx);
            status |= gr_poly_zero(T, ctx);
            status |= gr_poly_set_scalar(S, inv, ctx);
        }
        else if (lenB == 1)  /* lenA >= lenB = 1 */
        {
            status |= gr_inv(inv, B->coeffs, ctx);
            status |= gr_poly_set_scalar(T, inv, ctx);
            status |= gr_poly_one(G, ctx);
            status |= gr_poly_zero(S, ctx);
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

            status |= gr_inv(inv, GR_ENTRY(B->coeffs, lenB - 1, sz), ctx);
            status |= _gr_poly_xgcd_euclidean(&lenG, g, s, t, A->coeffs, lenA, B->coeffs, lenB, inv, ctx);

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

            if (gr_is_one(GR_ENTRY(G->coeffs, G->length - 1, sz), ctx) != T_TRUE)
            {
                status |= gr_inv(inv, GR_ENTRY(G->coeffs, G->length - 1, sz), ctx);
                status |= gr_poly_mul_scalar(G, G, inv, ctx);
                status |= gr_poly_mul_scalar(S, S, inv, ctx);
                status |= gr_poly_mul_scalar(T, T, inv, ctx);
            }
        }

        GR_TMP_CLEAR(inv, ctx);

        return status;
    }
}
