/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

/* todo: over appropriate domains, make_monic -> extract content; unit */

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

/* assumes lenA >= lenB >= 1, and both A and B have nonzero leading
   coefficient */
int
_gr_poly_gcd_euclidean(gr_ptr G, slong * lenG, gr_srcptr A, slong lenA,
                                gr_srcptr B, slong lenB, gr_ctx_t ctx)
{
    const slong lenW = FLINT_MAX(lenA - lenB + 1, lenB) + lenA + 2 * lenB;
    gr_ptr Q, R1, R2, R3, T, W;
    slong lenR2, lenR3;
    slong sz = ctx->sizeof_elem;
    int status;

    if (lenB == 1)
    {
        *lenG = 1;
        return gr_one(G, ctx);
    }

    status = GR_SUCCESS;

    /* lenA >= lenB > 1 */
    GR_TMP_INIT_VEC(W, lenW, ctx);
    Q = W;
    R1 = GR_ENTRY(W, FLINT_MAX(lenA - lenB + 1, lenB), sz);
    R2 = GR_ENTRY(R1, lenA, sz);
    R3 = GR_ENTRY(R2, lenB, sz);

    status |= _gr_poly_divrem(Q, R1, A, lenA, B, lenB, ctx);

    lenR3 = lenB - 1;
    status |= _gr_vec_normalise(&lenR3, R1, lenR3, ctx);

    if (status != GR_SUCCESS)
    {
        *lenG = 0;
        goto cleanup;
    }

    if (lenR3 == 0)
    {
        status |= _gr_vec_set(G, B, lenB, ctx);
        *lenG = lenB;
        goto cleanup;
    }

    T = R3;
    R3 = R1;
    R1 = T;
    status |= _gr_vec_set(R2, B, lenB, ctx);
    lenR2 = lenB;

    do
    {
        status |= _gr_poly_divrem(Q, R2, R2, lenR2, R3, lenR3, ctx);
        lenR2 = lenR3 - 1;
        status |= _gr_vec_normalise(&lenR2, R2, lenR2, ctx);

        if (status != GR_SUCCESS)
        {
            *lenG = 0;
            goto cleanup;
        }

        GR_VEC_SWAP(R2, lenR2, R3, lenR3);
    }
    while (lenR3 > 0);

    _gr_vec_swap(G, R2, lenR2, ctx);
    *lenG = lenR2;

cleanup:
    GR_TMP_CLEAR_VEC(W, lenW, ctx);

    return status;
}

int
gr_poly_gcd_euclidean(gr_poly_t G, const gr_poly_t A,
                        const gr_poly_t B, gr_ctx_t ctx)
{
    return gr_poly_gcd_wrapper((gr_method_poly_gcd_op) _gr_poly_gcd_euclidean, 1, G, A, B, ctx);
}

