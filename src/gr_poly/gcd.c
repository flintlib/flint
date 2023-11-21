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

#include "gr_poly.h"

/* todo: over appropriate domains, make_monic -> extract content; unit */
/* todo: extract powers of x (see fmpz_poly_gcd) */

/* assumes lenA >= lenB >= 1, and both A and B have nonzero leading
   coefficient */
int
_gr_poly_gcd(gr_ptr G, slong * lenG, gr_srcptr A, slong lenA,
                                gr_srcptr B, slong lenB, gr_ctx_t ctx)
{
    if (_gr_poly_gcd_euclidean(G, lenG, A, lenA, B, lenB, ctx) == GR_SUCCESS)
        return GR_SUCCESS;

    return GR_UNABLE;
}

int
gr_poly_gcd(gr_poly_t G, const gr_poly_t A,
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
        return gr_poly_gcd(G, B, A, ctx);

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

    status = _gr_poly_gcd(g, &lenG, A->coeffs, lenA, B->coeffs, lenB, ctx);

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
