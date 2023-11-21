/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

/* fixme: might return GR_DOMAIN where GR_UNABLE might be more
   appropriate, in cases where we have an unnormalized input
   over an inexact/uncomputable domain */

int
_gr_poly_make_monic(gr_ptr res, gr_srcptr poly, slong len, gr_ctx_t ctx)
{
    gr_srcptr lead;
    gr_ptr inv;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (len == 0)
        return GR_DOMAIN;

    lead = GR_ENTRY(poly, len - 1, sz);

    if (gr_is_one(lead, ctx) == T_TRUE)
    {
        status = _gr_vec_set(res, poly, len - 1, ctx);
    }
    else if (gr_is_neg_one(lead, ctx) == T_TRUE)
    {
        status = _gr_vec_neg(res, poly, len - 1, ctx);
    }
    else
    {
        GR_TMP_INIT(inv, ctx);

        /* try to multiply by inverse */
        status |= gr_inv(inv, lead, ctx);

        if (status == GR_SUCCESS)
        {
            status = _gr_vec_mul_scalar(res, poly, len - 1, inv, ctx);
        }
        else
        {
            /* otherwise try dividing */
            status = _gr_vec_div_scalar(res, poly, len - 1, lead, ctx);
        }

        GR_TMP_CLEAR(inv, ctx);
    }

    if (status == GR_SUCCESS)
        status = gr_one(GR_ENTRY(res, len - 1, sz), ctx);

    return status;
}

int
gr_poly_make_monic(gr_poly_t res, const gr_poly_t src, gr_ctx_t ctx)
{
    if (res != src)
    {
        int status;
        gr_poly_fit_length(res, src->length, ctx);
        status = _gr_poly_make_monic(res->coeffs, src->coeffs, src->length, ctx);
        _gr_poly_set_length(res, src->length, ctx);
        return status;
    }
    else
    {
        return _gr_poly_make_monic(res->coeffs, src->coeffs, src->length, ctx);
    }
}
