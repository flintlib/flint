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
acb_theta_agm_ext_step_bad(acb_ptr r, acb_srcptr a, acb_srcptr roots, slong g,
                           slong prec)
{
    slong k;

    for (k = 0; k < (1 << (g + 1)); k++)
    {
        acb_theta_agm_sqrt_lowprec(&r[k], &a[k], &roots[k], prec);
    }
    acb_theta_agm_ext_step_sqrt(r, r, g, prec);
}
