/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

void
_ca_poly_shift_right(ca_ptr res, ca_srcptr poly, slong len, slong n, ca_ctx_t ctx)
{
    slong i;

    /* Copy in forward order to avoid writing over unshifted coefficients */
    if (res != poly)
    {
        for (i = 0; i < len - n; i++)
            ca_set(res + i, poly + n + i, ctx);
    }
    else
    {
        for (i = 0; i < len - n; i++)
            ca_swap(res + i, res + n + i, ctx);
    }

}

void
ca_poly_shift_right(ca_poly_t res, const ca_poly_t poly, slong n, ca_ctx_t ctx)
{
    if (n == 0)
    {
        ca_poly_set(res, poly, ctx);
        return;
    }

    if (poly->length <= n)
    {
        ca_poly_zero(res, ctx);
        return;
    }

    ca_poly_fit_length(res, poly->length - n, ctx);
    _ca_poly_shift_right(res->coeffs, poly->coeffs, poly->length, n, ctx);
    _ca_poly_set_length(res, poly->length - n, ctx);
}
