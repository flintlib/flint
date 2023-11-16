/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

void
ca_poly_transfer(ca_poly_t res, ca_ctx_t res_ctx, const ca_poly_t src, ca_ctx_t src_ctx)
{
    slong i, len;

    if (res_ctx == src_ctx)
    {
        ca_poly_set(res, src, res_ctx);
    }
    else
    {
        len = src->length;

        ca_poly_fit_length(res, len, res_ctx);
        _ca_poly_set_length(res, len, res_ctx);

        for (i = 0; i < src->length; i++)
            ca_transfer(res->coeffs + i, res_ctx, src->coeffs + i, src_ctx);

        _ca_poly_normalise(res, res_ctx);
    }
}
