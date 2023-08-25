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
acb_theta_naive_ellipsoid(acb_theta_eld_t E, acb_ptr new_z, acb_ptr c, arb_ptr u,
    slong ord, acb_srcptr z, slong nb_z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong lp = ACB_THETA_LOW_PREC;
    arf_t R2;
    arf_t eps;
    arb_mat_t cho;
    arb_ptr offset;
    slong k;

    arf_init(R2);
    arf_init(eps);
    arb_mat_init(cho, g, g);
    offset = _arb_vec_init(g);

    acb_theta_eld_cho(cho, tau, prec);

    if (arb_mat_is_finite(cho))
    {
        /* Reduce all z, set offset and upper bounds */
        acb_theta_naive_radius(R2, eps, cho, ord, prec);
        acb_theta_naive_reduce(offset, new_z, c, u, z, nb_z, tau, cho, prec);
        for (k = 0; k < nb_z; k++)
        {
            arb_mul_arf(&u[k], &u[k], eps, prec);
        }

        /* Get radius for error of at most eps and fill ellipsoid */
        acb_theta_eld_fill(E, cho, R2, offset, lp);
    }
    else
    {
        /* Cannot compute cho, result will be nan */
        _acb_vec_zero(new_z, nb_z);
        arb_mat_one(cho);
        arf_zero(R2);
        acb_theta_eld_fill(E, cho, R2, offset, lp);
        for (k = 0; k < nb_z; k++)
        {
            acb_indeterminate(&c[k]);
            arb_pos_inf(&u[k]);
        }
    }

    arf_clear(R2);
    arf_clear(eps);
    arb_mat_clear(cho);
    _arb_vec_clear(offset, g);
}
