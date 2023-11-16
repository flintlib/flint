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
arb_atanh(arb_t z, const arb_t x, slong prec)
{
    if (arb_is_zero(x))
    {
        arb_zero(z);
    }
    else
    {
        arb_t t;
        arb_init(t);

        arb_sub_ui(t, x, 1, prec + 4);
        arb_div(t, x, t, prec + 4);
        arb_mul_2exp_si(t, t, 1);
        arb_neg(t, t);
        arb_log1p(z, t, prec);
        arb_mul_2exp_si(z, z, -1);

        arb_clear(t);
    }
}

