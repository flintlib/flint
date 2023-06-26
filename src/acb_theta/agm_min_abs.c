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
acb_theta_agm_min_abs(arb_t min, acb_srcptr a, slong nb, slong prec)
{
    arb_t abs;
    slong k;

    arb_init(abs);

    arb_pos_inf(min);
    for (k = 0; k < nb; k++)
    {
        acb_abs(abs, &a[k], prec);
        arb_min(min, min, abs, prec);
    }

    arb_clear(abs);
}
