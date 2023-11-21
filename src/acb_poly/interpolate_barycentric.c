/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_interpolate_barycentric(acb_ptr poly,
    acb_srcptr xs, acb_srcptr ys, slong n, slong prec)
{
    acb_ptr P, Q, w;
    acb_t t;
    slong i, j;

    if (n == 1)
    {
        acb_set(poly, ys);
        return;
    }

    P = _acb_vec_init(n + 1);
    Q = _acb_vec_init(n);
    w = _acb_vec_init(n);
    acb_init(t);

    _acb_poly_product_roots(P, xs, n, prec);

    for (i = 0; i < n; i++)
    {
        acb_one(w + i);

        for (j = 0; j < n; j++)
        {
            if (i != j)
            {
                acb_sub(t, xs + i, xs + j, prec);
                acb_mul(w + i, w + i, t, prec);
            }
        }

        acb_inv(w + i, w + i, prec);
    }

    _acb_vec_zero(poly, n);

    for (i = 0; i < n; i++)
    {
        _acb_poly_div_root(Q, t, P, n + 1, xs + i, prec);
        acb_mul(t, w + i, ys + i, prec);
        _acb_vec_scalar_addmul(poly, Q, n, t, prec);
    }

    _acb_vec_clear(P, n + 1);
    _acb_vec_clear(Q, n);
    _acb_vec_clear(w, n);
    acb_clear(t);
}

void
acb_poly_interpolate_barycentric(acb_poly_t poly,
    acb_srcptr xs, acb_srcptr ys, slong n, slong prec)
{
    if (n == 0)
    {
        acb_poly_zero(poly);
    }
    else
    {
        acb_poly_fit_length(poly, n);
        _acb_poly_set_length(poly, n);
        _acb_poly_interpolate_barycentric(poly->coeffs,
            xs, ys, n, prec);
        _acb_poly_normalise(poly);
    }
}
