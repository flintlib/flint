/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_theta_dupl_radius(arf_t rho, const arf_t r, acb_srcptr th, slong nb,
    slong prec)
{
    arb_t abs;
    arf_t bound, max;
    slong k;

    arb_init(abs);
    arf_init(bound);
    arf_init(max);

    arf_zero(max);
    for (k = 0; k < nb; k++)
    {
        acb_abs(abs, &th[k], prec);
        arb_get_ubound_arf(bound, abs, prec);
        arf_max(max, max, bound);
    }
    arf_div(rho, r, max, prec, ARF_RND_FLOOR);
    arf_mul_2exp_si(rho, rho, -1);
    arf_min(rho, rho, max);

    arb_clear(abs);
    arf_clear(bound);
    arf_clear(max);
}
