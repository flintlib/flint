/*
    Copyright (C) 2021 Matthias Gessinger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_graeffe_transform(acb_ptr b, acb_srcptr a, slong len, slong prec)
{
    slong lo, le, ls, deg, i;
    acb_ptr pe, po;

    if (len <= 1)
    {
        if (len)
            acb_sqr(b, a, prec);
        return;
    }

    deg = len - 1;
    lo = len / 2;
    ls = 2 * lo - 1;
    le = deg / 2 + 1;
    po = _acb_vec_init(lo);
    pe = _acb_vec_init(FLINT_MAX(le, ls));

    for (i = deg; i >= 0; i--)
    {
        if (i % 2 == 0)
            acb_set(pe + i / 2, a + i);
        else
            acb_set(po + i / 2, a + i);
    }

    _acb_poly_mul(b, pe, le, pe, le, prec);
    _acb_poly_mul(pe, po, lo, po, lo, prec);
    _acb_poly_sub(b + 1, b + 1, ls, pe, ls, prec);

    if (len % 2 == 0)
    {
        _acb_vec_neg(b, b, deg);
        acb_set(b + deg, pe + (deg - 1));
    }

    _acb_vec_clear(pe, FLINT_MAX(le, ls));
    _acb_vec_clear(po, lo);
}

void
acb_poly_graeffe_transform(acb_poly_t b, const acb_poly_t a, slong prec)
{
    acb_poly_fit_length(b, a->length);
    _acb_poly_graeffe_transform(b->coeffs, a->coeffs, a->length, prec);
    _acb_poly_set_length(b, a->length);
}
