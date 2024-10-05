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

#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"

#define GR_VEC_NORM(R, lenR) \
    (status) |= _gr_vec_normalise(&(lenR), (R), (lenR), (ctx))

#define __set(B, lenB, A, lenA)          \
do {                                     \
    status |= _gr_vec_set((B), (A), (lenA), ctx);     \
    (lenB) = (lenA);                     \
} while (0)

#define __add(C, lenC, A, lenA, B, lenB)                    \
do {                                                        \
    status |= _gr_poly_add((C), (A), (lenA), (B), (lenB), ctx); \
    (lenC) = FLINT_MAX((lenA), (lenB));                     \
    GR_VEC_NORM((C), (lenC));                             \
} while (0)

#define __sub(C, lenC, A, lenA, B, lenB)                    \
do {                                                        \
    status |= _gr_poly_sub((C), (A), (lenA), (B), (lenB), ctx); \
    (lenC) = FLINT_MAX((lenA), (lenB));                     \
    GR_VEC_NORM((C), (lenC));                             \
} while (0)

#define __mul(C, lenC, A, lenA, B, lenB)                            \
do {                                                                \
    if ((lenA) != 0 && (lenB) != 0)                                 \
    {                                                               \
        if ((lenA) >= (lenB))                                       \
            status |= _gr_poly_mul((C), (A), (lenA), (B), (lenB), ctx); \
        else                                                        \
            status |= _gr_poly_mul((C), (B), (lenB), (A), (lenA), ctx); \
        (lenC) = (lenA) + (lenB) - 1;                               \
    }                                                               \
    else                                                            \
    {                                                               \
        (lenC) = 0;                                                 \
    }                                                               \
} while (0)

#define __divrem(Q, lenQ, R, lenR, A, lenA, B, lenB)                          \
do {                                                                          \
    if ((lenA) >= (lenB))                                                     \
    {                                                                         \
        status |= _gr_poly_divrem((Q), (R), (A), (lenA), (B), (lenB), ctx);   \
        (lenQ) = (lenA) - (lenB) + 1;                                         \
        (lenR) = (lenB) - 1;                                                  \
        GR_VEC_NORM((R), (lenR));                                             \
    }                                                                         \
    else                                                                      \
    {                                                                         \
        status |= _gr_vec_set((R), (A), (lenA), ctx);                         \
        (lenQ) = 0;                                                           \
        (lenR) = (lenA);                                                      \
    }                                                                         \
} while (0)

/* todo: add _gr_poly_div */
#define __div(Q, lenQ, A, lenA, B, lenB)                                      \
do {                                                                          \
    if ((lenA) >= (lenB))                                                     \
    {                                                                         \
        status |= _gr_poly_div_newton((Q), (A), (lenA), (B), (lenB), ctx);    \
        (lenQ) = (lenA) - (lenB) + 1;                                         \
    }                                                                         \
    else                                                                      \
    {                                                                         \
        (lenQ) = 0;                                                           \
    }                                                                         \
} while (0)

int _gr_poly_xgcd_hgcd(slong * Glen, gr_ptr G, gr_ptr S, gr_ptr T,
                          gr_srcptr A, slong lenA, gr_srcptr B, slong lenB,
                          slong hgcd_cutoff,
                          slong cutoff,
                          gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong lenG, lenS, lenT;

    if (lenB == 1)
    {
        status |= gr_set(G, B, ctx);
        status |= gr_one(T, ctx);
        lenG = 1;
        lenS = 0;
        lenT = 1;
    }
    else
    {
        slong lenq, lenr, len1 = lenA + lenB;
        gr_ptr q, r, invB;

        GR_TMP_INIT_VEC(q, len1, ctx);
        r = GR_ENTRY(q, lenA, sz);
        GR_TMP_INIT(invB, ctx);

        __divrem(q, lenq, r, lenr, A, lenA, B, lenB);

        if (lenr == 0)
        {
            __set(G, lenG, B, lenB);
            status |= gr_one(T, ctx);
            lenS = 0;
            lenT = 1;
        }
        else
        {
            gr_ptr h, j, v, w, R[4], X;
            slong lenh, lenj, lenv, lenw, lenR[4], len2;
            slong sgnR;

            lenh = lenj = lenB;
            lenv = lenw = lenA + lenB - 2;
            lenR[0] = lenR[1] = lenR[2] = lenR[3] = (lenB + 1) / 2;

            len2 = 2 * lenh + 2 * lenv + 4 * lenR[0];

            GR_TMP_INIT_VEC(X, len2, ctx);
            h = X;
            j = GR_ENTRY(h, lenh, sz);
            v = GR_ENTRY(j, lenj, sz);
            w = GR_ENTRY(v, lenv, sz);
            R[0] = GR_ENTRY(w, lenw, sz);
            R[1] = GR_ENTRY(R[0], lenR[0], sz);
            R[2] = GR_ENTRY(R[1], lenR[1], sz);
            R[3] = GR_ENTRY(R[2], lenR[2], sz);

            status |= _gr_poly_hgcd(NULL, &sgnR, R, lenR, h, &lenh, j, &lenj, B, lenB, r, lenr, hgcd_cutoff, ctx);

            if (sgnR > 0)
            {
                status |= _gr_vec_neg(S, R[1], lenR[1], ctx);
                status |= _gr_vec_set(T, R[0], lenR[0], ctx);
            }
            else
            {
                status |= _gr_vec_set(S, R[1], lenR[1], ctx);
                status |= _gr_vec_neg(T, R[0], lenR[0], ctx);
            }
            lenS = lenR[1];
            lenT = lenR[0];

            while (lenj != 0 && status == GR_SUCCESS)
            {
                __divrem(q, lenq, r, lenr, h, lenh, j, lenj);
                __mul(v, lenv, q, lenq, T, lenT);

                {
                    slong l;
                    _gr_vec_swap(S, T, FLINT_MAX(lenS, lenT), ctx);
                    l = lenS; lenS = lenT; lenT = l;
                }
                if (lenr != 0) /* prevent overflow of T on last iteration */
                   __sub(T, lenT, T, lenT, v, lenv);
                else /* lenr == 0 */
                {
                    __set(G, lenG, j, lenj);

                    goto cofactor;
                }

                if (lenj < cutoff)
                {
                    gr_ptr u0 = R[0], u1 = R[1];
                    slong lenu0 = lenr - 1, lenu1 = lenj - 1;

                    status |= _gr_poly_xgcd_euclidean(&lenG, G, u0, u1, j, lenj, r, lenr, ctx);

                    GR_VEC_NORM(u0, lenu0);
                    GR_VEC_NORM(u1, lenu1);

                    __mul(v, lenv, S, lenS, u0, lenu0);
                    __mul(w, lenw, T, lenT, u1, lenu1);
                    __add(S, lenS, v, lenv, w, lenw);

                    goto cofactor;
                }

                status |= _gr_poly_hgcd(NULL, &sgnR, R, lenR, h, &lenh, j, &lenj, j,lenj, r, lenr, hgcd_cutoff, ctx);

                __mul(v, lenv, R[1], lenR[1], T, lenT);
                __mul(w, lenw, R[2], lenR[2], S, lenS);

                __mul(q, lenq, S, lenS, R[3], lenR[3]);
                if (sgnR > 0)
                    __sub(S, lenS, q, lenq, v, lenv);
                else
                    __sub(S, lenS, v, lenv, q, lenq);

                __mul(q, lenq, T, lenT, R[0], lenR[0]);
                if (sgnR > 0)
                    __sub(T, lenT, q, lenq, w, lenw);
                else
                    __sub(T, lenT, w, lenw, q, lenq);
            }
            __set(G, lenG, h, lenh);

            cofactor:

            __mul(v, lenv, S, lenS, A, lenA);
            __sub(w, lenw, G, lenG, v, lenv);
            __div(T, lenT, w, lenw, B, lenB);

            GR_TMP_CLEAR_VEC(X, len2, ctx);
        }

        GR_TMP_CLEAR_VEC(q, len1, ctx);
        GR_TMP_CLEAR(invB, ctx);
    }

    status |= _gr_vec_zero(GR_ENTRY(S, lenS, sz), lenB - 1 - lenS, ctx);
    status |= _gr_vec_zero(GR_ENTRY(T, lenT, sz), lenA - 1 - lenT, ctx);

    *Glen = lenG;

    return status;
}

int
gr_poly_xgcd_hgcd(gr_poly_t G, gr_poly_t S, gr_poly_t T, const gr_poly_t A, const gr_poly_t B, slong hgcd_cutoff, slong cutoff, gr_ctx_t ctx)
{
    if (A->length < B->length)
    {
        return gr_poly_xgcd_hgcd(G, T, S, B, A, hgcd_cutoff, cutoff, ctx);
    }
    else                        /* lenA >= lenB >= 0 */
    {
        const slong lenA = A->length, lenB = B->length;
        int status = GR_SUCCESS;
        slong sz = ctx->sizeof_elem;

        if (lenA == 0)          /* lenA = lenB = 0 */
        {
            status |= gr_poly_zero(G, ctx);
            status |= gr_poly_zero(S, ctx);
            status |= gr_poly_zero(T, ctx);
        }
        else if (lenB == 0)     /* lenA > lenB = 0 */
        {
            gr_ptr t;

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
            gr_ptr t;

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

            status |= _gr_poly_xgcd_hgcd(&lenG, g, s, t, A->coeffs, lenA, B->coeffs, lenB, hgcd_cutoff, cutoff, ctx);
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
