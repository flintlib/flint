/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

#define CA_VEC_NORM(success, R, lenR, ctx) \
    do { \
        (success) = 1; \
        while ((lenR) > 0) \
        { \
            truth_t is_zero; \
            is_zero = ca_check_is_zero((R) + (lenR) - 1, (ctx)); \
            if (is_zero == T_TRUE) \
                (lenR)--; \
            else if (is_zero == T_UNKNOWN) \
            { \
                (success) = 0; \
                break; \
            } \
            else \
            { \
                break; \
            } \
        } \
    } while (0) \

/* assumes lenA >= lenB >= 1, and both A and B have nonzero leading
   coefficient */
slong
_ca_poly_gcd_euclidean(ca_ptr G, ca_srcptr A, slong lenA,
                                ca_srcptr B, slong lenB, ca_ctx_t ctx)
{
    const slong lenW = FLINT_MAX(lenA - lenB + 1, lenB) + lenA + 2 * lenB;
    ca_t invR3;
    ca_ptr Q, R1, R2, R3, T, W;
    slong lenR2, lenR3;
    int success;

    if (lenB == 1)
    {
        ca_one(G, ctx);
        return 1;
    }

    success = 1;

    /* lenA >= lenB > 1 */
    ca_init(invR3, ctx);
    W = _ca_vec_init(lenW, ctx);
    Q = W;
    R1 = W + FLINT_MAX(lenA - lenB + 1, lenB);
    R2 = R1 + lenA;
    R3 = R2 + lenB;

    ca_inv(invR3, B + lenB - 1, ctx);
    _ca_poly_divrem(Q, R1, A, lenA, B, lenB, invR3, ctx);

    lenR3 = lenB - 1;
    CA_VEC_NORM(success, R1, lenR3, ctx);
    if (!success)
        goto cleanup;

    if (lenR3 == 0)
    {
        ca_clear(invR3, ctx);
        _ca_vec_set(G, B, lenB, ctx);
        _ca_vec_clear(W, lenW, ctx);
        return lenB;
    }

    T = R3;
    R3 = R1;
    R1 = T;
    _ca_vec_set(R2, B, lenB, ctx);
    lenR2 = lenB;

    do
    {
        ca_inv(invR3, R3 + (lenR3 - 1), ctx);
        _ca_poly_divrem(Q, R1, R2, lenR2, R3, lenR3, invR3, ctx);

        lenR2 = lenR3--;
        CA_VEC_NORM(success, R1, lenR3, ctx);
        if (!success)
            goto cleanup;

        T = R2;
        R2 = R3;
        R3 = R1;
        R1 = T;
    }
    while (lenR3 > 0);

    _ca_vec_set(G, R2, lenR2, ctx);

cleanup:
    _ca_vec_clear(W, lenW, ctx);
    ca_clear(invR3, ctx);

    if (success)
        return lenR2;
    else
        return 0;
}

int
ca_poly_gcd_euclidean(ca_poly_t G, const ca_poly_t A,
                        const ca_poly_t B, ca_ctx_t ctx)
{
    slong lenA = A->length, lenB = B->length, lenG;
    ca_ptr g;

    if (A->length == 0 && B->length == 0)
    {
        ca_poly_zero(G, ctx);
        return 1;
    }

    if (A->length == 0)
        return ca_poly_make_monic(G, B, ctx);

    if (B->length == 0)
        return ca_poly_make_monic(G, A, ctx);

    if (A->length < B->length)
        return ca_poly_gcd_euclidean(G, B, A, ctx);

    if (ca_check_is_zero(A->coeffs + A->length - 1, ctx) != T_FALSE ||
        ca_check_is_zero(B->coeffs + B->length - 1, ctx) != T_FALSE)
    {
        return 0;
    }

    /* lenA >= lenB >= 1 */
    if (G == A || G == B)
    {
        g = _ca_vec_init(FLINT_MIN(lenA, lenB), ctx);
    }
    else
    {
        ca_poly_fit_length(G, FLINT_MIN(lenA, lenB), ctx);
        g = G->coeffs;
    }

    lenG = _ca_poly_gcd_euclidean(g, A->coeffs, lenA, B->coeffs, lenB, ctx);

    if (G == A || G == B)
    {
        _ca_vec_clear(G->coeffs, G->alloc, ctx);
        G->coeffs = g;
        G->alloc = FLINT_MIN(lenA, lenB);
        G->length = FLINT_MIN(lenA, lenB);
    }
    _ca_poly_set_length(G, lenG, ctx);

    if (lenG == 0)
    {
        return 0;
    }
    else
    {
        if (G->length == 1)
            ca_one(G->coeffs, ctx);
        else
            ca_poly_make_monic(G, G, ctx);
        return 1;
    }
}
