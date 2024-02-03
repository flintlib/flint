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
#include "gr_mat.h"

int _gr_poly_resultant_sylvester(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx)
{
    slong d, e, i;
    int status = GR_SUCCESS;
    gr_mat_t M;

    d = len1 - 1;
    e = len2 - 1;

    gr_mat_init(M, d + e, d + e, ctx);

    for (i = 0; i < e; i++)
        status |= _gr_poly_reverse(gr_mat_entry_ptr(M, i, i, ctx), poly1, len1, len1, ctx);

    for (i = 0; i < d; i++)
        status |= _gr_poly_reverse(gr_mat_entry_ptr(M, e + i, i, ctx), poly2, len2, len2, ctx);

    status |= gr_mat_det(res, M, ctx);

    gr_mat_clear(M, ctx);

    return status;
}

int gr_poly_resultant_sylvester(gr_ptr r, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx)
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
        status |= _gr_poly_resultant_sylvester(r, f->coeffs, len1,  g->coeffs, len2, ctx);
    }
    else
    {
        status |= _gr_poly_resultant_sylvester(r, g->coeffs, len2, f->coeffs, len1, ctx);

        if (((len1 | len2) & 1) == 0)
            status |= gr_neg(r, r, ctx);
    }

    return status;
}
