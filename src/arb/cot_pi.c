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
arb_cot_pi(arb_t y, const arb_t x, slong prec)
{

    if (!arb_is_finite(x))
    {
        arb_indeterminate(y);
    }
    else if (arb_is_int_2exp_si(x, -1))
    {
        if (arb_is_int(x))
            arb_indeterminate(y);
        else
            arb_zero(y);
    }
    else
    {
        arb_t u;
        arb_init(u);
        arb_sin_cos_pi(y, u, x, prec + 4);
        arb_div(y, u, y, prec);
        arb_clear(u);
    }
}

