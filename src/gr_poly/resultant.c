/*
    Copyright (C) 2007, 2008 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

/* Important: fmpz_mod_poly currently relies on these tuning values.
   If the are changed to accommodate other rings, fmpz_mod_poly_resultant
   should override the tuning values. */
#define HGCD_CUTOFF 200
#define HGCD_INNER_CUTOFF 100

int _gr_poly_resultant(gr_ptr res, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (_gr_poly_resultant_small(res, A, lenA, B, lenB, ctx) == GR_SUCCESS)
        return GR_SUCCESS;

    if (FLINT_MIN(lenA, lenB) >= HGCD_CUTOFF && gr_ctx_is_finite(ctx) == T_TRUE)
        status = _gr_poly_resultant_hgcd(res, A, lenA, B, lenB, HGCD_INNER_CUTOFF, HGCD_CUTOFF, ctx);
    else
        status = _gr_poly_resultant_euclidean(res, A, lenA, B, lenB, ctx);

    /* A division-free algorithm should succeed */
    if (status != GR_SUCCESS)
        status = _gr_poly_resultant_sylvester(res, A, lenA, B, lenB, ctx);

    return status;
}

int gr_poly_resultant(gr_ptr r, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx)
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
        status |= _gr_poly_resultant(r, f->coeffs, len1,  g->coeffs, len2, ctx);
    }
    else
    {
        status |= _gr_poly_resultant(r, g->coeffs, len2, f->coeffs, len1, ctx);

        if (((len1 | len2) & 1) == 0)
            status |= gr_neg(r, r, ctx);
    }

    return status;
}
