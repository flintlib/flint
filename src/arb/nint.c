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
arb_nint(arb_t res, const arb_t x, slong prec)
{
    if (arb_is_int(x))
    {
        arb_set(res, x);
    }
    else
    {
        arb_t t, u;
        arb_init(t);
        arb_init(u);

        arb_set_d(t, 0.5);
        arb_add(t, x, t, prec);

        arb_mul_2exp_si(u, x, 1);
        arb_sub_ui(u, u, 1, prec);
        arb_mul_2exp_si(u, u, -2);

        arb_floor(res, t, prec);

        /* nint(x) = floor(x+0.5) - isint((2*x-1)/4) */

        if (arb_is_int(u))
        {
            arb_sub_ui(res, res, 1, prec);
        }
        else if (arb_contains_int(u))
        {
            arf_one(arb_midref(u));
            mag_one(arb_radref(u));
            arb_mul_2exp_si(u, u, -1);
            arb_sub_ui(res, res, 1, prec);
        }
    }
}
