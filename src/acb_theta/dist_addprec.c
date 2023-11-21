/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "acb_theta.h"

slong
acb_theta_dist_addprec(const arb_t d2)
{
    arb_t x;
    slong prec = ACB_THETA_LOW_PREC;
    slong res;

    arb_init(x);
    arb_const_log2(x, prec);
    arb_div(x, d2, x, prec);

    if (arb_is_finite(x) && (arf_cmpabs_2exp_si(arb_midref(x), 30) <= 0))
    {
        res = arf_get_si(arb_midref(x), prec);
    }
    else /* should never happen */
    {
        res = 0;
    }

    arb_clear(x);
    return res;
}
