/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "acb_theta.h"

slong
acb_theta_sum_addprec(const arb_t d2)
{
    arb_t x;
    arf_t b;
    slong prec = ACB_THETA_LOW_PREC;
    slong res;

    arb_init(x);
    arf_init(b);

    arb_const_log2(x, prec);
    arb_div(x, d2, x, prec);
    arb_get_ubound_arf(b, x, prec);

    if (arf_is_finite(b) && (arf_cmpabs_2exp_si(b, 40) <= 0))
    {
        res = arf_get_si(b, prec);
    }
    else /* should never happen */
    {
        res = 0;
    }

    arb_clear(x);
    return res;
}
