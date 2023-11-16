/*
    Copyright (C) 2021 Matthias Gessinger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_graeffe_transform(arb_ptr b, arb_srcptr a, slong len, slong prec)
{
    slong lo, le, ls, deg, i;
    arb_ptr pe, po;

    if (len <= 1)
    {
        if (len)
            arb_sqr(b, a, prec);
        return;
    }

    deg = len - 1;
    lo = len / 2;
    ls = 2 * lo - 1;
    le = deg / 2 + 1;
    po = _arb_vec_init(lo);
    pe = _arb_vec_init(FLINT_MAX(le, ls));

    for (i = deg; i >= 0; i--)
    {
        if (i % 2 == 0)
            arb_set(pe + i / 2, a + i);
        else
            arb_set(po + i / 2, a + i);
    }

    _arb_poly_mul(b, pe, le, pe, le, prec);
    _arb_poly_mul(pe, po, lo, po, lo, prec);
    _arb_poly_sub(b + 1, b + 1, ls, pe, ls, prec);

    if (len % 2 == 0)
    {
        _arb_vec_neg(b, b, deg);
        arb_set(b + deg, pe + (deg - 1));
    }

    _arb_vec_clear(pe, FLINT_MAX(le, ls));
    _arb_vec_clear(po, lo);
}

void
arb_poly_graeffe_transform(arb_poly_t b, const arb_poly_t a, slong prec)
{
    arb_poly_fit_length(b, a->length);
    _arb_poly_graeffe_transform(b->coeffs, a->coeffs, a->length, prec);
    _arb_poly_set_length(b, a->length);
}
