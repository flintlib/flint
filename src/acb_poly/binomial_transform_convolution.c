/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "acb_poly.h"

void
_acb_poly_binomial_transform_convolution(acb_ptr b, acb_srcptr a, slong alen, slong len, slong prec)
{
    slong i;
    acb_ptr c, d;

    alen = FLINT_MIN(alen, len);

    c = _acb_vec_init(alen);
    d = _acb_vec_init(len);

    _acb_poly_borel_transform(c, a, alen, prec);
    for (i = 1; i < alen; i += 2)
        acb_neg(c + i, c + i);

    acb_one(d);
    for (i = 1; i < len; i++)
        acb_div_ui(d + i, d + i - 1, i, prec);

    _acb_poly_mullow(b, d, len, c, alen, len, prec);

    _acb_poly_inv_borel_transform(b, b, len, prec);

    _acb_vec_clear(c, alen);
    _acb_vec_clear(d, len);
}

void
acb_poly_binomial_transform_convolution(acb_poly_t b, const acb_poly_t a, slong len, slong prec)
{
    if (len == 0 || a->length == 0)
    {
        acb_poly_zero(b);
        return;
    }

    if (b == a)
    {
        acb_poly_t c;
        acb_poly_init2(c, len);
        _acb_poly_binomial_transform_convolution(c->coeffs, a->coeffs, a->length, len, prec);
        acb_poly_swap(b, c);
        acb_poly_clear(c);
    }
    else
    {
        acb_poly_fit_length(b, len);
        _acb_poly_binomial_transform_convolution(b->coeffs, a->coeffs, a->length, len, prec);
    }

    _acb_poly_set_length(b, len);
    _acb_poly_normalise(b);
}

