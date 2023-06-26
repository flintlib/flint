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
acb_theta_agm_rel_dist(arb_t eps, acb_srcptr a, slong nb, slong lowprec, slong prec)
{
    arb_t abs;

    arb_init(abs);

    acb_theta_agm_abs_dist(eps, a, nb, lowprec, prec);
    acb_abs(abs, &a[0], lowprec);
    arb_div(eps, eps, abs, lowprec);

    arb_clear(abs);
}
