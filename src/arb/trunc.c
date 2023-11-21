/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_trunc(arb_t res, const arb_t x, slong prec)
{
    if (arb_contains_zero(x))
    {
        arb_t a;
        arb_init(a);

        mag_one(arb_radref(a));

        if (arb_contains(a, x))
        {
            arb_zero(res);
        }
        else
        {
            arb_t b;
            arb_init(b);
            arb_floor(a, x, prec);
            arb_ceil(b, x, prec);
            arb_union(res, a, b, prec);
            arb_clear(b);
        }

        arb_clear(a);
    }
    else if (arf_sgn(arb_midref(x)) > 0)
    {
        arb_floor(res, x, prec);
    }
    else
    {
        arb_ceil(res, x, prec);
    }
}
