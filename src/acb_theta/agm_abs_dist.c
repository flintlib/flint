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
acb_theta_agm_abs_dist(arb_t eps, acb_srcptr a, slong nb, slong lowprec,
                       slong prec)
{
    acb_t diff;
    arb_t abs;
    slong k;

    acb_init(diff);
    arb_init(abs);

    arb_zero(eps);
    for (k = 1; k < nb; k++)
    {
        acb_sub(diff, &a[k], &a[0], prec);
        acb_abs(abs, diff, lowprec);
        arb_max(eps, eps, abs, lowprec);
    }

    acb_clear(diff);
    arb_clear(abs);
}
