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
_acb_poly_taylor_shift_divconquer(acb_ptr poly, const acb_t c, slong len, slong prec)
{
    acb_struct d[2];

    if (len <= 1 || acb_is_zero(c))
        return;

    if (len == 2)
    {
        acb_addmul(poly, poly + 1, c, prec);
        return;
    }

    d[0] = *c;

    acb_init(d + 1);
    acb_one(d + 1); /* no need to free */

    _acb_poly_compose_divconquer(poly, poly, len, d, 2, prec);
}

void
acb_poly_taylor_shift_divconquer(acb_poly_t g, const acb_poly_t f,
    const acb_t c, slong prec)
{
    if (f != g)
        acb_poly_set_round(g, f, prec);

    _acb_poly_taylor_shift_divconquer(g->coeffs, c, g->length, prec);
}

