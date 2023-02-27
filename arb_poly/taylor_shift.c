/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_taylor_shift(arb_ptr poly, const arb_t c, slong n, slong prec)
{
    if (n <= 30 || (n <= 500 && arb_bits(c) == 1 && n < 30 + 3 * sqrt(prec))
                || (n <= 100 && arb_bits(c) < 0.01 * prec))
    {
        _arb_poly_taylor_shift_horner(poly, c, n, prec);
    }
    else if (prec > 2 * n)
    {
        _arb_poly_taylor_shift_convolution(poly, c, n, prec);
    }
    else
    {
        _arb_poly_taylor_shift_divconquer(poly, c, n, prec);
    }
}

void
arb_poly_taylor_shift(arb_poly_t g, const arb_poly_t f,
    const arb_t c, slong prec)
{
    if (f != g)
        arb_poly_set_round(g, f, prec);

    _arb_poly_taylor_shift(g->coeffs, c, g->length, prec);
}

