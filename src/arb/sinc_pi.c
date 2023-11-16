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
arb_sinc_pi(arb_t res, const arb_t x, slong prec)
{
    mag_t m;
    arb_t t;

    if (!arb_is_finite(x))
    {
        if (arf_is_nan(arb_midref(x)))
            arb_indeterminate(res);
        else if (arb_contains_zero(x))
            arb_zero_pm_one(res);
        else
            arb_zero(res);
        return;
    }

    if (arb_is_int(x))
    {
        if (arb_is_zero(x))
            arb_one(res);
        else
            arb_zero(res);
        return;
    }

    mag_init(m);
    arb_init(t);

    arb_get_mag_lower(m, x);

    if (mag_cmp_2exp_si(m, -1) > 0)
    {
        arb_const_pi(t, prec + 4);
        arb_mul(t, t, x, prec + 4);
        arb_sin_pi(res, x, prec + 4);
        arb_div(res, res, t, prec);
    }
    else
    {
        arb_const_pi(t, prec + 4);
        arb_mul(t, t, x, prec + 4);
        arb_sinc(res, t, prec);
    }

    mag_clear(m);
    arb_clear(t);
}

