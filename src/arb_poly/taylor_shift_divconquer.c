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
_arb_poly_taylor_shift_divconquer(arb_ptr poly, const arb_t c, slong len, slong prec)
{
    arb_struct d[2];

    if (len <= 1 || arb_is_zero(c))
        return;

    if (len == 2)
    {
        arb_addmul(poly, poly + 1, c, prec);
        return;
    }

    d[0] = *c;

    arb_init(d + 1);
    arb_one(d + 1); /* no need to free */

    _arb_poly_compose_divconquer(poly, poly, len, d, 2, prec);
}

void
arb_poly_taylor_shift_divconquer(arb_poly_t g, const arb_poly_t f,
    const arb_t c, slong prec)
{
    if (f != g)
        arb_poly_set_round(g, f, prec);

    _arb_poly_taylor_shift_divconquer(g->coeffs, c, g->length, prec);
}

