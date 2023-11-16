/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_binomial_transform(acb_ptr b, acb_srcptr a, slong alen, slong len, slong prec)
{
    if (alen < 10 || len < 10)
        _acb_poly_binomial_transform_basecase(b, a, alen, len, prec);
    else
        _acb_poly_binomial_transform_convolution(b, a, alen, len, prec);
}

void
acb_poly_binomial_transform(acb_poly_t b, const acb_poly_t a, slong len, slong prec)
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
        _acb_poly_binomial_transform(c->coeffs, a->coeffs, a->length, len, prec);
        acb_poly_swap(b, c);
        acb_poly_clear(c);
    }
    else
    {
        acb_poly_fit_length(b, len);
        _acb_poly_binomial_transform(b->coeffs, a->coeffs, a->length, len, prec);
    }

    _acb_poly_set_length(b, len);
    _acb_poly_normalise(b);
}

