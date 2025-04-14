/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

/* todo: defer to make_monic or explicitly 1 the leading coefficient
         when over a field */
int
_gr_poly_canonicalise_unit(gr_ptr res, gr_srcptr poly, slong len, gr_ctx_t ctx)
{
    gr_srcptr lead;
    gr_ptr inv;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (len == 0)
        return GR_SUCCESS;

    if (gr_ctx_is_field(ctx) == T_TRUE)
        return _gr_poly_make_monic(res, poly, len, ctx);

    lead = GR_ENTRY(poly, len - 1, sz);

    if (gr_is_one(lead, ctx) == T_TRUE)
    {
        if (res != poly)
            status = _gr_vec_set(res, poly, len, ctx);
    }
    else if (gr_is_neg_one(lead, ctx) == T_TRUE)
    {
        status = _gr_vec_neg(res, poly, len, ctx);
    }
    else if (gr_is_zero(lead, ctx) != T_FALSE)
    {
        return GR_UNABLE;
    }
    else
    {
        GR_TMP_INIT(inv, ctx);

        status |= gr_canonical_unit(inv, lead, ctx);
        status |= gr_inv(inv, inv, ctx);

        if (status == GR_SUCCESS)
        {
            status = _gr_vec_mul_scalar(res, poly, len, inv, ctx);
        }

        GR_TMP_CLEAR(inv, ctx);
    }

    return status;
}

int
gr_poly_canonicalise_unit(gr_poly_t res, const gr_poly_t src, gr_ctx_t ctx)
{
    if (res != src)
    {
        int status;
        gr_poly_fit_length(res, src->length, ctx);
        status = _gr_poly_canonicalise_unit(res->coeffs, src->coeffs, src->length, ctx);
        _gr_poly_set_length(res, src->length, ctx);
        _gr_poly_normalise(res, ctx);
        return status;
    }
    else
    {
        return _gr_poly_canonicalise_unit(res->coeffs, src->coeffs, src->length, ctx);
    }
}

