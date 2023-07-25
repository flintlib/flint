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
acb_theta_ql_step_1(acb_ptr r, acb_srcptr th0, acb_srcptr th, acb_srcptr roots,
    arb_srcptr dist0, arb_srcptr dist, slong g, slong prec)
{
    slong n = 1 << g;

    acb_theta_agm_mul_tight(r, th0, th, dist0, dist, g, prec);
    _acb_vec_scalar_mul_2exp_si(r, r, n, g);
    acb_theta_agm_sqrt(r, r, roots, n, prec);
}
