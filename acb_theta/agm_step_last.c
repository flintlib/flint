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
acb_theta_agm_step_last(acb_t r, acb_srcptr a, slong g, slong prec)
{
    slong k;
    slong n = 1 << g;

    acb_zero(r);
    for (k = 0; k < n; k++)
    {
        acb_add(r, r, &a[k], prec);
    }
    acb_mul_2exp_si(r, r, -g);
}
