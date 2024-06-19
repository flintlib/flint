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

#define GR_VEC_NORM(status, R, lenR, sz, ctx) \
    do { \
        (void) sz; \
        (status) |= _gr_vec_normalise(&(lenR), (R), (lenR), (ctx)); \
    } while (0)

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

int
_gr_poly_xgcd_euclidean(slong * lenG, gr_ptr G, gr_ptr S, gr_ptr T, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
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

        /* todo: preinv1 B; use also for second division later */
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
                status |= _gr_poly_divrem(Q, D, D, lenD, V3, lenV3, ctx);
                lenQ = lenD - lenV3 + 1;
                lenD = lenV3 - 1;
                GR_VEC_NORM(status, D, lenD, sz, ctx);
                if (status != GR_SUCCESS)
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
                GR_VEC_SWAP(D, lenD, V3, lenV3);

            } while (lenV3 != 0);

            status |= _gr_vec_set(G, D, lenD, ctx);
            status |= _gr_vec_set(S, U, lenU, ctx);

            if (status == GR_SUCCESS)
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
    return gr_poly_xgcd_wrapper((gr_method_poly_xgcd_op) _gr_poly_xgcd_euclidean, G, S, T, A, B, ctx);
}
