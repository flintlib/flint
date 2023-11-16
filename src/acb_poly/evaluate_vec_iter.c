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
_acb_poly_evaluate_vec_iter(acb_ptr ys, acb_srcptr poly, slong plen,
    acb_srcptr xs, slong n, slong prec)
{
    slong i;

    for (i = 0; i < n; i++)
        _acb_poly_evaluate(ys + i, poly, plen, xs + i, prec);
}

void
acb_poly_evaluate_vec_iter(acb_ptr ys,
        const acb_poly_t poly, acb_srcptr xs, slong n, slong prec)
{
    _acb_poly_evaluate_vec_iter(ys, poly->coeffs,
                                        poly->length, xs, n, prec);
}
