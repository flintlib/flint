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
acb_theta_g2_chi6m2(acb_poly_t res, acb_srcptr dth, slong prec)
{
    slong g = 2;
    slong n = 1 << (2 * g);
    slong nb = acb_theta_jet_nb(1, g);
    acb_ptr th;
    acb_t den;
    slong k;

    th = _acb_vec_init(n);
    acb_init(den);

    for (k = 0; k < n; k++)
    {
        acb_set(&th[k], &dth[k * nb]);
    }
    acb_theta_g2_chi63(res, dth, prec);
    acb_theta_g2_chi5(den, th, prec);
    acb_poly_scalar_div(res, res, den, prec);

    _acb_vec_clear(th, n);
    acb_clear(den);
}
