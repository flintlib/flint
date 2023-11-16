/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_binomial_transform(arb_ptr b, arb_srcptr a, slong alen, slong len, slong prec)
{
    if (alen < 10 || len < 10)
        _arb_poly_binomial_transform_basecase(b, a, alen, len, prec);
    else
        _arb_poly_binomial_transform_convolution(b, a, alen, len, prec);
}

void
arb_poly_binomial_transform(arb_poly_t b, const arb_poly_t a, slong len, slong prec)
{
    if (len == 0 || a->length == 0)
    {
        arb_poly_zero(b);
        return;
    }

    if (b == a)
    {
        arb_poly_t c;
        arb_poly_init2(c, len);
        _arb_poly_binomial_transform(c->coeffs, a->coeffs, a->length, len, prec);
        arb_poly_swap(b, c);
        arb_poly_clear(c);
    }
    else
    {
        arb_poly_fit_length(b, len);
        _arb_poly_binomial_transform(b->coeffs, a->coeffs, a->length, len, prec);
    }

    _arb_poly_set_length(b, len);
    _arb_poly_normalise(b);
}

