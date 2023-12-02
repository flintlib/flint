/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"
#include "gr_vec.h"

#define GR_VEC_NORM(status, R, lenR, sz, ctx) \
    do { \
        (void) sz; \
        (status) |= _gr_vec_normalise(&(lenR), (R), (lenR), (ctx)); \
    } while (0)


#define __set(B, lenB, A, lenA)                     \
do {                                                \
    status |= _gr_vec_set((B), (A), (lenA), ctx);   \
    (lenB) = (lenA);                                \
} while (0)

#define __divrem(Q, lenQ, R, lenR, A, lenA, B, lenB)                \
do {                                                                \
    if ((lenA) >= (lenB))                                           \
    {                                                               \
        status |= _gr_poly_divrem((Q), (R), (A), (lenA), (B), (lenB), ctx); \
        (lenQ) = (lenA) - (lenB) + 1;                               \
        (lenR) = (lenB) - 1;                                        \
        GR_VEC_NORM(status, R, lenR, sz, ctx);                      \
    }                                                               \
    else                                                            \
    {                                                               \
        status |= _gr_vec_set((R), (A), (lenA), ctx);               \
        (lenQ) = 0;                                                 \
        (lenR) = (lenA);                                            \
    }                                                               \
} while (0)

/*
    XXX: Incidentally, this implementation currently supports aliasing.
    But since this may change in the future, no function other than
    gr_poly_gcd_hgcd() should rely on this.
 */

int
_gr_poly_gcd_hgcd(gr_ptr G, slong * _lenG, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong inner_cutoff, slong cutoff, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong lenG, lenJ, lenR, lenQ;
    slong sz = ctx->sizeof_elem;
    slong alloc = 2 * lenB + lenA;
    gr_ptr J, R, Q;

    if (lenB > lenA)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    GR_TMP_INIT_VEC(J, alloc, ctx);
    R = GR_ENTRY(J, lenB, sz);
    Q = GR_ENTRY(R, lenB, sz);

    __divrem(Q, lenQ, R, lenR, A, lenA, B, lenB);
    (void) lenQ; /* unused */

    if (lenR == 0)
    {
        __set(G, lenG, B, lenB);
    }
    else
    {
        status |= _gr_poly_hgcd(NULL, NULL, NULL, NULL, G, &(lenG), J, &(lenJ), B, lenB, R, lenR, inner_cutoff, ctx);

        while (lenJ != 0)
        {
            /* todo: only rem */
            __divrem(Q, lenQ, R, lenR, G, lenG, J, lenJ);
            (void) lenQ; /* unused */

            if (lenR == 0)
            {
                __set(G, lenG, J, lenJ);
                break;
            }

            if (lenJ < cutoff)
            {
                status |= _gr_poly_gcd_euclidean(G, &lenG, J, lenJ, R, lenR, ctx);
                break;
            }

            status |= _gr_poly_hgcd(NULL, NULL, NULL, NULL, G, &(lenG), J, &(lenJ), J, lenJ, R, lenR, inner_cutoff, ctx);
        }
    }

    GR_TMP_CLEAR_VEC(J, alloc, ctx);
    *_lenG = lenG;

    return status;
}

int gr_poly_gcd_hgcd(gr_poly_t G, const gr_poly_t A, const gr_poly_t B, slong inner_cutoff, slong cutoff, gr_ctx_t ctx)
{
    if (A->length < B->length)
    {
        return gr_poly_gcd_hgcd(G, B, A, inner_cutoff, cutoff, ctx);
    }
    else                        /* lenA >= lenB >= 0 */
    {
        slong lenA = A->length, lenB = B->length, lenG;
        gr_poly_t tG;
        gr_ptr g;
        int status = GR_SUCCESS;

        if (lenA == 0)          /* lenA = lenB = 0 */
        {
            return gr_poly_zero(G, ctx);
        }
        else if (lenB == 0)     /* lenA > lenB = 0 */
        {
            return gr_poly_make_monic(G, A, ctx);
        }
        else                    /* lenA >= lenB >= 1 */
        {
            if (G == A || G == B)
            {
                gr_poly_init2(tG, FLINT_MIN(lenA, lenB), ctx);
                g = tG->coeffs;
            }
            else
            {
                gr_poly_fit_length(G, FLINT_MIN(lenA, lenB), ctx);
                g = G->coeffs;
            }

            status = _gr_poly_gcd_hgcd(g, &lenG, A->coeffs, lenA, B->coeffs, lenB, inner_cutoff, cutoff, ctx);

            if (G == A || G == B)
            {
                gr_poly_swap(tG, G, ctx);
                gr_poly_clear(tG, ctx);
            }
            G->length = lenG;

            if (G->length == 1)
                status |= gr_one(G->coeffs, ctx);
            else
                status |= gr_poly_make_monic(G, G, ctx);

            return status;
        }
    }
}
