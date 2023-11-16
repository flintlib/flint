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
arb_coth(arb_t y, const arb_t x, slong prec)
{
    arb_t t, u;
    int sign = arf_sgn(arb_midref(x)) < 0;

    arb_init(t);
    arb_init(u);

    arb_mul_2exp_si(t, x, 1);

    if (!sign)
        arb_neg(t, t);

    if (arf_cmpabs_2exp_si(arb_midref(x), 1) > 0)
    {
        /* coth(x) = 1 + 2 exp(-2x) / (1 - exp(-2x)) */
        arb_exp(t, t, prec + 4);
        arb_sub_ui(u, t, 1, prec + 4);
        arb_div(y, t, u, prec + 4);
        arb_mul_2exp_si(y, y, 1);
        arb_sub_ui(y, y, 1, prec);
    }
    else
    {
        /* coth(x) = (exp(2x) + 1) / (exp(2x) - 1) */
        arb_expm1(t, t, prec + 4);
        arb_add_ui(y, t, 2, prec + 4);
        arb_div(y, y, t, prec);
    }

    if (!sign)
        arb_neg(y, y);

    arb_clear(t);
    arb_clear(u);
}
