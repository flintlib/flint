/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

/* assumes lenA >= lenB >= 1, and both A and B have nonzero leading
   coefficient */
slong
_ca_poly_gcd(ca_ptr G, ca_srcptr A, slong lenA,
                                ca_srcptr B, slong lenB, ca_ctx_t ctx)
{
    return _ca_poly_gcd_euclidean(G, A, lenA, B, lenB, ctx);
}

int
ca_poly_gcd(ca_poly_t G, const ca_poly_t A,
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
        return ca_poly_gcd(G, B, A, ctx);

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

    lenG = _ca_poly_gcd(g, A->coeffs, lenA, B->coeffs, lenB, ctx);

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
