/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_acosh(arb_t z, const arb_t x, slong prec)
{
    if (arb_is_one(x))
    {
        arb_zero(z);
    }
    else
    {
        arb_t t;
        arb_init(t);

        arb_mul(t, x, x, prec + 4);
        arb_sub_ui(t, t, 1, prec + 4);
        arb_sqrt(t, t, prec + 4);
        arb_add(t, t, x, prec + 4);
        arb_log(z, t, prec);

        arb_clear(t);
    }
}

