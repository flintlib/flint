/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_csch(arb_t res, const arb_t x, slong prec)
{
    if (arb_contains_zero(x) || !arb_is_finite(x))
    {
        arb_indeterminate(res);
    }
    else if (arf_cmpabs_2exp_si(arb_midref(x), 0) > 0)
    {
        arb_t t;
        arb_init(t);

        if (arf_sgn(arb_midref(x)) > 0)
        {
            arb_neg(t, x);
            arb_exp(t, t, prec + 4);
            arb_mul(res, t, t, prec + 4);
            arb_sub_ui(res, res, 1, prec + 4);
            arb_div(res, t, res, prec);
            arb_neg(res, res);
        }
        else
        {
            arb_exp(t, x, prec + 4);
            arb_mul(res, t, t, prec + 4);
            arb_sub_ui(res, res, 1, prec + 4);
            arb_div(res, t, res, prec);
        }

        arb_mul_2exp_si(res, res, 1);
        arb_clear(t);
    }
    else
    {
        arb_sinh(res, x, prec + 4);
        arb_inv(res, res, prec);
    }
}

