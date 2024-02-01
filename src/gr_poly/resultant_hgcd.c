/*
    Copyright (C) 2007, 2008 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

#define GR_VEC_NORM(status, R, lenR, sz, ctx) \
    do { \
        (void) sz; \
        (status) |= _gr_vec_normalise(&(lenR), (R), (lenR), (ctx)); \
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
    gr_poly_resultant_hgcd() should rely on this.
 */

int _gr_poly_resultant_hgcd(gr_ptr res, gr_srcptr A, slong lenA,
                               gr_srcptr B, slong lenB, slong inner_cutoff, slong cutoff, gr_ctx_t ctx)
{
    slong len1 = FLINT_MIN(lenA, lenB), len2 = 2 * lenB;
    gr_ptr G, J, R, Q, lc;
    slong lenG, lenJ, lenR, lenQ;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    GR_TMP_INIT_VEC(G, len1, ctx);
    GR_TMP_INIT_VEC(J, len2, ctx);
    GR_TMP_INIT(lc, ctx);

    R = GR_ENTRY(J, lenB, sz);

    /* todo: avoid Q alloc when we have rem */
    GR_TMP_INIT_VEC(Q, lenA, ctx);

    status |= gr_one(res, ctx);
    status |= gr_set(lc, GR_ENTRY(B, lenB - 1, sz), ctx);

    /* todo: only rem */
    __divrem(Q, lenQ, R, lenR, A, lenA, B, lenB);
    (void) lenQ;

    if (status != GR_SUCCESS)
        goto cleanup;

    if (lenR == 0)
    {
        if (lenB == 1)
        {
            status |= gr_pow_ui(lc, lc, lenA - 1, ctx);
            status |= gr_mul(res, res, lc, ctx);
        }
        else
            status |= gr_zero(res, ctx);
    }
    else
    {
        status |= gr_pow_ui(lc, lc, lenA - lenR, ctx);
        status |= gr_mul(res, res, lc, ctx);

        if (((lenA | lenB) & 1) == 0)
            status |= gr_neg(res, res, ctx);

        status |= _gr_poly_hgcd(res, NULL, NULL, NULL, G, &(lenG), J, &(lenJ), B, lenB, R, lenR, inner_cutoff, ctx);

        if (status != GR_SUCCESS)
            goto cleanup;

        while (lenJ != 0)
        {
            status |= gr_set(lc, GR_ENTRY(J, lenJ - 1, sz), ctx);

            /* todo: only rem */
            __divrem(Q, lenQ, R, lenR, G, lenG, J, lenJ);
            (void) lenQ;

            if (status != GR_SUCCESS)
                goto cleanup;

            if (lenR == 0)
            {
                if (lenJ == 1)
                {
                    status |= gr_pow_ui(lc, lc, lenG - 1, ctx);
                    status |= gr_mul(res, res, lc, ctx);
                }
                else
                {
                    status |= gr_zero(res, ctx);
                }

                break;
            }
            else
            {
                status |= gr_pow_ui(lc, lc, lenG - lenR, ctx);
                status |= gr_mul(res, res, lc, ctx);

                if (((lenG | lenJ) & 1) == 0)
                    status |= gr_neg(res, res, ctx);
            }

            if (lenJ < cutoff)
            {
                gr_ptr r;
                GR_TMP_INIT(r, ctx);
                status |= _gr_poly_resultant_euclidean(r, J, lenJ, R, lenR, ctx);
                status |= gr_mul(res, res, r, ctx);
                GR_TMP_CLEAR(r, ctx);
                break;
            }

            status |= _gr_poly_hgcd(res, NULL, NULL, NULL, G, &(lenG), J, &(lenJ), J, lenJ, R, lenR, inner_cutoff, ctx);
        }
    }

cleanup:
    GR_TMP_CLEAR_VEC(J, len2, ctx);
    GR_TMP_CLEAR_VEC(G, len1, ctx);
    GR_TMP_CLEAR(lc, ctx);

    GR_TMP_CLEAR_VEC(Q, lenA, ctx);

    return status;
}

int gr_poly_resultant_hgcd(gr_ptr r, const gr_poly_t f,
                             const gr_poly_t g, slong inner_cutoff, slong cutoff, gr_ctx_t ctx)
{
    slong len1 = f->length;
    slong len2 = g->length;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (len1 == 0 || len2 == 0)
    {
        return gr_zero(r, ctx);
    }

    if (gr_is_zero(GR_ENTRY(f->coeffs, len1 - 1, sz), ctx) != T_FALSE ||
        gr_is_zero(GR_ENTRY(g->coeffs, len2 - 1, sz), ctx) != T_FALSE)
    {
        return GR_UNABLE;
    }

    if (len1 >= len2)
    {
        status |= _gr_poly_resultant_hgcd(r, f->coeffs, len1,  g->coeffs, len2, inner_cutoff, cutoff, ctx);
    }
    else
    {
        status |= _gr_poly_resultant_hgcd(r, g->coeffs, len2, f->coeffs, len1, inner_cutoff, cutoff, ctx);

        if (((len1 | len2) & 1) == 0)
            status |= gr_neg(r, r, ctx);
    }

    return status;
}
