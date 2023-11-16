/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_product_roots(acb_ptr poly, acb_srcptr xs, slong n, slong prec)
{
    if (n == 0)
    {
        acb_one(poly);
    }
    else if (n == 1)
    {
        acb_neg(poly, xs);
        acb_one(poly + 1);
    }
    else if (n == 2)
    {
        acb_mul(poly, xs + 0, xs + 1, prec);
        acb_add(poly + 1, xs + 0, xs + 1, prec);
        acb_neg(poly + 1, poly + 1);
        acb_one(poly + 2);
    }
    else if (n == 3)
    {
        acb_mul(poly + 1, xs, xs + 1, prec);
        acb_mul(poly, poly + 1, xs + 2, prec);
        acb_neg(poly, poly);
        acb_add(poly + 2, xs, xs + 1, prec);
        acb_addmul(poly + 1, poly + 2, xs + 2, prec);
        acb_add(poly + 2, poly + 2, xs + 2, prec);
        acb_neg(poly + 2, poly + 2);
        acb_one(poly + 3);
    }
    else
    {
        const slong m = (n + 1) / 2;
        acb_ptr tmp;

        tmp = _acb_vec_init(n + 2);

        _acb_poly_product_roots(tmp, xs, m, prec);
        _acb_poly_product_roots(tmp + m + 1, xs + m, n - m, prec);
        _acb_poly_mul_monic(poly, tmp, m + 1, tmp + m + 1, n - m + 1, prec);

        _acb_vec_clear(tmp, n + 2);
    }
}

void
acb_poly_product_roots(acb_poly_t poly, acb_srcptr xs, slong n, slong prec)
{
    acb_poly_fit_length(poly, n + 1);
    _acb_poly_product_roots(poly->coeffs, xs, n, prec);
    _acb_poly_set_length(poly, n + 1);
}
