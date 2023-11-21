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
_ca_poly_set_roots(ca_ptr poly, ca_srcptr xs, const ulong * exp, slong n, ca_ctx_t ctx)
{
    if (n == 0 || (n == 1 && exp[0] == 0))
    {
        ca_one(poly, ctx);
    }
    else if (n == 1)
    {
        if (exp[0] == 1)
        {
            ca_neg(poly, xs, ctx);
            ca_one(poly + 1, ctx);
        }
        else if (exp[0] == 2)
        {
            ca_sqr(poly, xs, ctx);
            ca_mul_si(poly + 1, xs, -2, ctx);
            ca_one(poly + 2, ctx);
        }
        else
        {
            slong i, e;

            e = exp[0];

            ca_one(poly + e, ctx);
            for (i = e - 1; i >= 0; i--)
            {
                ca_mul(poly + i, poly + i + 1, xs, ctx);
                ca_mul_si(poly + i, poly + i, -(i + 1), ctx);
                ca_div_ui(poly + i, poly + i, e - i, ctx);
            }
        }
    }
    else if (n == 2 && exp[0] == 1 && exp[1] == 1)
    {
        ca_mul(poly, xs + 0, xs + 1, ctx);
        ca_add(poly + 1, xs + 0, xs + 1, ctx);
        ca_neg(poly + 1, poly + 1, ctx);
        ca_one(poly + 2, ctx);
    }
    else
    {
        slong i, m, deg_left, deg_right;
        ca_ptr tmp;

        /* todo: balance degrees */
        m = (n + 1) / 2;

        deg_left = deg_right = 0;
        for (i = 0; i < m; i++)
            deg_left += exp[i];
        for (i = m; i < n; i++)
            deg_right += exp[i];

        tmp = _ca_vec_init(deg_left + deg_right + 2, ctx);

        _ca_poly_set_roots(tmp, xs, exp, m, ctx);
        _ca_poly_set_roots(tmp + deg_left + 1, xs + m, exp + m, n - m, ctx);
        _ca_poly_mul(poly, tmp, deg_left + 1, tmp + deg_left + 1, deg_right + 1, ctx);

        _ca_vec_clear(tmp, deg_left + deg_right + 2, ctx);
    }
}

void
ca_poly_set_roots(ca_poly_t poly, ca_vec_t roots, const ulong * exp, ca_ctx_t ctx)
{
    slong i, deg, len;

    len = ca_vec_length(roots, ctx);

    /* todo: check for overflow? */
    deg = 0;
    for (i = 0; i < len; i++)
        deg += exp[i];

    ca_poly_fit_length(poly, deg + 1, ctx);
    _ca_poly_set_roots(poly->coeffs, roots->entries, exp, len, ctx);
    _ca_poly_set_length(poly, deg + 1, ctx);
}
