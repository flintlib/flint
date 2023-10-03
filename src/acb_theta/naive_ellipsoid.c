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
acb_theta_naive_ellipsoid(acb_theta_eld_t E, acb_ptr new_zs, acb_ptr cs, arb_ptr us,
    acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arf_t R2;
    arf_t eps;
    arb_mat_t C;
    arb_ptr v;
    slong k;

    arf_init(R2);
    arf_init(eps);
    arb_mat_init(C, g, g);
    v = _arb_vec_init(g);

    acb_theta_eld_cho(C, tau, prec);

    if (arb_mat_is_finite(C))
    {
        /* Get radius, fill ellipsoid */
        acb_theta_naive_radius(R2, eps, C, 0, prec);
        acb_theta_naive_reduce(v, new_zs, cs, us, zs, nb, tau, C, prec);
        for (k = 0; k < nb; k++)
        {
            arb_mul_arf(&us[k], &us[k], eps, prec);
        }
        acb_theta_eld_fill(E, C, R2, v);
    }
    else
    {
        /* Cannot compute C, result will be nan */
        _acb_vec_zero(new_zs, nb);
        arb_mat_one(C);
        arf_zero(R2);
        acb_theta_eld_fill(E, C, R2, v);
        for (k = 0; k < nb; k++)
        {
            acb_indeterminate(&cs[k]);
            arb_pos_inf(&us[k]);
        }
    }

    arf_clear(R2);
    arf_clear(eps);
    arb_mat_clear(C);
    _arb_vec_clear(v, g);
}
