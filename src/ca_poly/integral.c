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
_ca_poly_integral(ca_ptr res, ca_srcptr poly, slong len, ca_ctx_t ctx)
{
    slong k = len - 1;

    for (k = len - 1; k > 0; k--)
        ca_div_ui(res + k, poly + k - 1, k, ctx);

    ca_zero(res, ctx);
}

void
ca_poly_integral(ca_poly_t res, const ca_poly_t poly, ca_ctx_t ctx)
{
    ca_poly_fit_length(res, poly->length + 1, ctx);
    _ca_poly_integral(res->coeffs, poly->coeffs, poly->length + 1, ctx);
    _ca_poly_set_length(res, poly->length + 1, ctx);
    _ca_poly_normalise(res, ctx);
}
