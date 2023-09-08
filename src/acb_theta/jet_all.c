/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

static void
acb_theta_jet_all_mid(acb_ptr dth, slong ord, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong b = ord + 1;
    slong hprec;
    slong lp = ACB_THETA_LOW_PREC;
    slong nb_jet = acb_theta_jet_nb(ord, g + 1);
    arb_t c, rho, t;
    arf_t eps, err;
    acb_ptr zetas, new_z, all_val, val, jet;
    slong k, kmod, j;

    arf_init(eps);
    arf_init(err);
    arb_init(c);
    arb_init(rho);
    arb_init(t);
    zetas = _acb_vec_init(b);
    new_z = _acb_vec_init(g);
    all_val = _acb_vec_init(n * n * n_pow(b, g));
    val = _acb_vec_init(n_pow(b, g));
    jet = _acb_vec_init(nb_jet);

    /* Get bounds and precision */
    acb_theta_jet_bounds_1(c, rho, z, tau, ord, lp);
    acb_theta_jet_radius(eps, err, c, rho, ord, g, prec, lp);
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
        _acb_vec_add(new_z, new_z, z, g, hprec);
        acb_theta_all(all_val + k * n * n, new_z, tau, 0, hprec);
    }

    /* Call jet_fd on each theta_{a,b} */
    for (k = 0; k < n * n; k++)
    {
        for (j = 0; j < n_pow(b, g); j++)
        {
            acb_set(&val[j], &all_val[j * n * n + k]);
        }
        acb_theta_jet_fd(jet, eps, err, val, ord, g, prec);
        for (j = 0; j < nb_jet; j++)
        {
            acb_set(&dth[j * n * n + k], &jet[j]);
        }
    }

    arf_clear(eps);
    arf_clear(err);
    arb_clear(c);
    arb_clear(rho);
    arb_clear(t);
    _acb_vec_clear(zetas, b);
    _acb_vec_clear(new_z, g);
    _acb_vec_clear(all_val, n * n * n_pow(b, g));
    _acb_vec_clear(val, n_pow(b, g));
    _acb_vec_clear(jet, nb_jet);
}

void
acb_theta_jet_all(acb_ptr dth, slong ord, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong nb = acb_theta_jet_nb(ord, g + 1);
    slong lp = ACB_THETA_LOW_PREC;
    acb_mat_t w;
    acb_ptr x;
    arb_t c, rho, b;
    arf_t eps, t;
    fmpz_t m;
    slong k, j;

    acb_mat_init(w, g, g);
    x = _acb_vec_init(g);
    arb_init(c);
    arb_init(rho);
    arb_init(b);
    arf_init(eps);
    arf_init(t);
    fmpz_init(m);

    /* Get bounds */
    acb_theta_jet_bounds_2(c, rho, z, tau, prec);

    /* Call jet_all_mid on midpoint */
    acb_mat_get_mid(w, tau);
    for (k = 0; k < g; k++)
    {
        acb_get_mid(&x[k], &z[k]);
    }
    acb_theta_jet_all_mid(dth, ord, x, w, prec);

    /* Add error bounds */
    arf_zero(eps);
    for (k = 0; k < g; k++)
    {
        arf_set_mag(t, arb_radref(acb_realref(&z[k])));
        arf_max(eps, eps, t);
        arf_set_mag(t, arb_radref(acb_imagref(&z[k])));
        arf_max(eps, eps, t);
        for (j = 0; j < g; j++)
        {
            arf_set_mag(t, arb_radref(acb_realref(acb_mat_entry(tau, k, j))));
            arf_max(eps, eps, t);
            arf_set_mag(t, arb_radref(acb_imagref(acb_mat_entry(tau, k, j))));
            arf_max(eps, eps, t);
        }
    }
    arb_get_lbound_arf(t, rho, lp);

    if (arf_cmp(t, eps) <= 0)
    {
        _acb_vec_indeterminate(dth, nb);
    }
    else
    {
        arb_sub_arf(b, rho, eps, lp);
        arb_pow_ui(b, b, ord + 1, lp);
        arb_div(b, c, b, lp);
        fmpz_fac_ui(m, ord + 1);
        arb_mul_fmpz(b, b, m, lp);
        fmpz_bin_uiui(m, ord + 1 + g + (g * g + g)/2, ord + 1);
        arb_mul_fmpz(b, b, m, lp);
        arb_mul_arf(b, b, eps, prec);
        for (k = 0; k < nb; k++)
        {
            acb_add_error_arb(&dth[k], b);
        }
    }

    acb_mat_clear(w);
    _acb_vec_clear(x, g);
    arb_clear(c);
    arb_clear(rho);
    arb_clear(b);
    arf_clear(eps);
    arf_clear(t);
    fmpz_clear(m);

}
