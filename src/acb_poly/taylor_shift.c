/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_taylor_shift(acb_ptr poly, const acb_t c, slong n, slong prec)
{
    if (n <= 30 || (n <= 500 && acb_bits(c) == 1 && n < 30 + 3 * sqrt(prec))
                || (n <= 100 && acb_bits(c) < 0.01 * prec))
    {
        _acb_poly_taylor_shift_horner(poly, c, n, prec);
    }
    else if (prec > 2 * n)
    {
        _acb_poly_taylor_shift_convolution(poly, c, n, prec);
    }
    else
    {
        _acb_poly_taylor_shift_divconquer(poly, c, n, prec);
    }
}

void
acb_poly_taylor_shift(acb_poly_t g, const acb_poly_t f,
    const acb_t c, slong prec)
{
    if (f != g)
        acb_poly_set_round(g, f, prec);

    _acb_poly_taylor_shift(g->coeffs, c, g->length, prec);
}

