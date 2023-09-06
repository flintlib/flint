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
acb_theta_jet_all(acb_ptr dth, slong ord, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong b = ord + 1;
    slong hprec;
    slong lp = ACB_THETA_LOW_PREC;
    slong nb_jet = acb_theta_jet_nb(ord, g + 1);
    arb_t c, rho, t;
    arf_t eps;
    acb_ptr zetas, new_z, all_val, val, jet;
    slong k, kmod, j;

    arf_init(eps);
    arb_init(c);
    arb_init(rho);
    arb_init(t);
    zetas = _acb_vec_init(b);
    new_z = _acb_vec_init(g);
    all_val = _acb_vec_init(n * n * n_pow(b, g));
    val = _acb_vec_init(n_pow(b, g));
    jet = _acb_vec_init(nb_jet);

    /* Get bounds and precision */
    acb_theta_jet_bounds(c, rho, z, tau, ord, lp);
    acb_theta_jet_radius(eps, c, rho, ord, g, prec, lp);
    arb_set_arf(t, eps);
    arb_log_base_ui(t, t, 2, lp);
    arb_neg(t, t);
    hprec = prec + arf_get_si(arb_midref(t), ARF_RND_CEIL);

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
        arb_set_arf(t, eps);
        _acb_vec_scalar_mul_arb(new_z, new_z, g, t, hprec);
        _acb_vec_add(new_z, new_z, z, g, prec); /* todo: need to get mid and adjust errors */
        acb_theta_all(all_val + k * n * n, new_z, tau, 0, prec);
    }

    /* Call jet_fd on each theta_{a,b} */
    for (k = 0; k < n * n; k++)
    {
        for (j = 0; j < n_pow(b, g); j++)
        {
            acb_set(&val[j], &all_val[j * n * n + k]);
        }
        acb_theta_jet_fd(jet, eps, c, rho, val, ord, g, prec);
        for (j = 0; j < nb_jet; j++)
        {
            acb_set(&dth[j * n * n + k], &jet[j]);
        }
    }

    arf_clear(eps);
    arb_clear(c);
    arb_clear(rho);
    arb_clear(t);
    _acb_vec_clear(zetas, b);
    _acb_vec_clear(new_z, g);
    _acb_vec_clear(all_val, n * n * n_pow(b, g));
    _acb_vec_clear(val, n_pow(b, g));
    _acb_vec_clear(jet, nb_jet);
}
