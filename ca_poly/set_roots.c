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
_ca_poly_set_roots(ca_ptr poly, ca_srcptr xs, slong n, ca_ctx_t ctx)
{
    if (n == 0)
    {
        ca_one(poly, ctx);
    }
    else if (n == 1)
    {
        ca_neg(poly, xs, ctx);
        ca_one(poly + 1, ctx);
    }
    else if (n == 2)
    {
        ca_mul(poly, xs + 0, xs + 1, ctx);
        ca_add(poly + 1, xs + 0, xs + 1, ctx);
        ca_neg(poly + 1, poly + 1, ctx);
        ca_one(poly + 2, ctx);
    }
    else
    {
        const slong m = (n + 1) / 2;
        ca_ptr tmp;

        tmp = _ca_vec_init(n + 2, ctx);

        _ca_poly_set_roots(tmp, xs, m, ctx);
        _ca_poly_set_roots(tmp + m + 1, xs + m, n - m, ctx);
        _ca_poly_mul(poly, tmp, m + 1, tmp + m + 1, n - m + 1, ctx);

        _ca_vec_clear(tmp, n + 2, ctx);
    }
}

void
ca_poly_set_roots(ca_poly_t poly, ca_vec_t roots, ca_ctx_t ctx)
{
    slong len = ca_vec_length(roots, ctx);
    ca_poly_fit_length(poly, len + 1, ctx);
    _ca_poly_set_roots(poly->coeffs, roots->coeffs, len, ctx);
    _ca_poly_set_length(poly, len + 1, ctx);
}

