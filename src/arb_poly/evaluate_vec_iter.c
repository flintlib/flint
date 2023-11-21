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
_arb_poly_evaluate_vec_iter(arb_ptr ys, arb_srcptr poly, slong plen,
    arb_srcptr xs, slong n, slong prec)
{
    slong i;

    for (i = 0; i < n; i++)
        _arb_poly_evaluate(ys + i, poly, plen, xs + i, prec);
}

void
arb_poly_evaluate_vec_iter(arb_ptr ys,
        const arb_poly_t poly, arb_srcptr xs, slong n, slong prec)
{
    _arb_poly_evaluate_vec_iter(ys, poly->coeffs,
                                        poly->length, xs, n, prec);
}
