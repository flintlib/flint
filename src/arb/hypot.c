/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_hypot(arb_t z, const arb_t x, const arb_t y, slong prec)
{
    if (arb_is_zero(y))
    {
        arb_abs(z, x);
    }
    else if (arb_is_zero(x))
    {
        arb_abs(z, y);
    }
    else
    {
        arb_t t;
        arb_init(t);
        /* TODO: use arb_fmma */
        arb_mul(t, x, x, prec + 4);
        arb_mul(z, y, y, prec + 4);
        arb_add(t, t, z, prec + 4);
        arb_sqrtpos(z, t, prec);
        arb_clear(t);
    }
}

