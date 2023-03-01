/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

void
ca_poly_set_fmpq_poly(ca_poly_t res, const fmpq_poly_t src, ca_ctx_t ctx)
{
    slong i;
    ca_poly_fit_length(res, src->length, ctx);
    if (fmpz_is_one(fmpq_poly_denref(src)))
    {
        for (i = 0; i < src->length; i++)
            ca_set_fmpz(res->coeffs + i, src->coeffs + i, ctx);
    }
    else
    {
        for (i = 0; i < src->length; i++)
        {
            ca_set_fmpz(res->coeffs + i, src->coeffs + i, ctx);
            ca_div_fmpz(res->coeffs + i, res->coeffs + i, fmpq_poly_denref(src), ctx);
        }
    }
    _ca_poly_set_length(res, src->length, ctx);
}
