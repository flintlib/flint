/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"
#include "gr_vec.h"

int
gr_poly_truncate(gr_poly_t poly, const gr_poly_t src, slong newlen, gr_ctx_t ctx)
{
    if (src == poly)
    {
        if (poly->length > newlen)
        {
            _gr_poly_set_length(poly, newlen, ctx);
            _gr_poly_normalise(poly, ctx);
        }

        return GR_SUCCESS;
    }
    else
    {
        int status = GR_SUCCESS;
        slong len = src->length;

        newlen = FLINT_MIN(newlen, src->length);

        gr_poly_fit_length(poly, newlen, ctx);
        status |= _gr_vec_set(poly->coeffs, src->coeffs, newlen, ctx);
        _gr_poly_set_length(poly, newlen, ctx);
        if (newlen < len)
            _gr_poly_normalise(poly, ctx);
        return status;
    }
}
