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
acb_theta_0b_jet(acb_ptr dth, slong ord, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong b = ord + 1;
    slong hprec;
    slong lp = ACB_THETA_LOW_PREC;
    arb_t eps, c, rho, t;
    acb_ptr zetas, new_z, all_val, val, all_jets;
    slong k, kmod, j;

    arb_init(eps);
    arb_init(c);
    arb_init(rho);
    arb_init(t);
    zetas = _acb_vec_init(b);
    new_z = _acb_vec_init(g);
    all_val = _acb_vec_init(n * n_pow(b, g));
    val = _acb_vec_init(n_pow(b, g));
    all_jets = _acb_vec_init(n * acb_theta_deriv_nb(ord, g + 1));

    /* Get bounds and precision */
    acb_theta_deriv_bounds(eps, c, rho, z, tau, ord, prec, lp);
    arb_log_base_ui(t, eps, 2, lp);
    arb_neg(t, t);
    hprec = prec + arf_get_si(t, ARF_RND_CEIL);

    /* Get all values */
    _acb_vec_unit_roots(zetas, b, b, hprec);
    for (k = 0; k < n_pow(b, g); k++)
    {
        kmod = k;
        for (j = 0; j < g; j++)
        {
            acb_set(&new_z[j], &zetas[kmod % b]);
            kmod = kmod / b;
        }
        _acb_vec_scalar_mul_arb(new_z, new_z, g, eps, hprec);
        _acb_vec_add(new_z, new_z, z, g, prec); /* todo: need to get mid and adjust errors */
        acb_theta_
    }
    
}
