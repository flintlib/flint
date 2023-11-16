/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_interpolate_barycentric(arb_ptr poly,
    arb_srcptr xs, arb_srcptr ys, slong n, slong prec)
{
    arb_ptr P, Q, w;
    arb_t t;
    slong i, j;

    if (n == 1)
    {
        arb_set(poly, ys);
        return;
    }

    P = _arb_vec_init(n + 1);
    Q = _arb_vec_init(n);
    w = _arb_vec_init(n);
    arb_init(t);

    _arb_poly_product_roots(P, xs, n, prec);

    for (i = 0; i < n; i++)
    {
        arb_one(w + i);

        for (j = 0; j < n; j++)
        {
            if (i != j)
            {
                arb_sub(t, xs + i, xs + j, prec);
                arb_mul(w + i, w + i, t, prec);
            }
        }

        arb_inv(w + i, w + i, prec);
    }

    _arb_vec_zero(poly, n);

    for (i = 0; i < n; i++)
    {
        _arb_poly_div_root(Q, t, P, n + 1, xs + i, prec);
        arb_mul(t, w + i, ys + i, prec);
        _arb_vec_scalar_addmul(poly, Q, n, t, prec);
    }

    _arb_vec_clear(P, n + 1);
    _arb_vec_clear(Q, n);
    _arb_vec_clear(w, n);
    arb_clear(t);
}

void
arb_poly_interpolate_barycentric(arb_poly_t poly,
    arb_srcptr xs, arb_srcptr ys, slong n, slong prec)
{
    if (n == 0)
    {
        arb_poly_zero(poly);
    }
    else
    {
        arb_poly_fit_length(poly, n);
        _arb_poly_set_length(poly, n);
        _arb_poly_interpolate_barycentric(poly->coeffs,
            xs, ys, n, prec);
        _arb_poly_normalise(poly);
    }
}
