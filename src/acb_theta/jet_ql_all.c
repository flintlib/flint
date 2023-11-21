/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "acb_mat.h"
#include "acb_theta.h"

static void
acb_theta_jet_ql_all_red(acb_ptr dth, acb_srcptr z, const acb_mat_t tau, slong ord, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n2 = 1 << (2 * g);
    slong b = ord + 1;
    slong hprec;
    slong lp = ACB_THETA_LOW_PREC;
    slong nb = acb_theta_jet_nb(ord, g);
    slong nb_low = acb_theta_jet_nb(ord + 2, g);
    int hasz = !_acb_vec_is_zero(z, g);
    arb_t c, rho, t;
    arf_t eps, err, e;
    acb_mat_t tau_mid;
    acb_ptr z_mid, zetas, new_z, all_val, val, jet, dth_low;
    arb_ptr err_vec;
    slong k, kmod, j;

    arb_init(c);
    arb_init(rho);
    arb_init(t);
    arf_init(eps);
    arf_init(err);

    /* Get bounds and high precision, fail if too large */
    acb_theta_jet_ql_bounds(c, rho, z, tau, ord);
    acb_theta_jet_ql_radius(eps, err, c, rho, ord, g, prec);
    arb_set_arf(t, eps);
    arb_log_base_ui(t, t, 2, lp);
    arb_neg(t, t);

    /* we expect that the second bound holds if finite, but still check */
    if (!arb_is_finite(t) || (arf_cmpabs_2exp_si(arb_midref(t), 20) > 0))
    {
        _acb_vec_indeterminate(dth, n2 * nb);
        arb_clear(c);
        arb_clear(rho);
        arb_clear(t);
        arf_clear(eps);
        arf_clear(err);
        return;
    }

    arf_init(e);
    acb_mat_init(tau_mid, g, g);
    z_mid = _acb_vec_init(g);
    zetas = _acb_vec_init(b);
    new_z = _acb_vec_init(g);
    all_val = _acb_vec_init(n2 * n_pow(b, g));
    val = _acb_vec_init(n_pow(b, g));
    jet = _acb_vec_init(nb);
    dth_low = _acb_vec_init(n2 * nb_low);
    err_vec = _arb_vec_init(nb);

    hprec = prec + ord * (arf_get_si(arb_midref(t), ARF_RND_CEIL) + g);
    arf_one(e);
    arf_mul_2exp_si(e, e, -hprec);

    /* Get midpoint at precision hprec */
    for (j = 0; j < g; j++)
    {
        for (k = 0; k < g; k++)
        {
            acb_get_mid(acb_mat_entry(tau_mid, j, k), acb_mat_entry(tau, j, k));
            acb_add_error_arf(acb_mat_entry(tau_mid, j, k), e);
            acb_set(acb_mat_entry(tau_mid, k, j), acb_mat_entry(tau_mid, j, k));
        }
        acb_get_mid(&z_mid[j], &z[j]);
        if (hasz)
        {
            acb_add_error_arf(&z_mid[j], e);
        }
    }

    /* Collect values around midpoint */
    _acb_vec_unit_roots(zetas, b, b, hprec);
    for (k = 0; k < n_pow(b, g); k++)
    {
        kmod = k;
        for (j = g - 1; j >= 0; j--)
        {
            acb_set(&new_z[j], &zetas[kmod % b]);
            kmod = kmod / b;
        }
        arb_set_arf(t, eps);
        _acb_vec_scalar_mul_arb(new_z, new_z, g, t, hprec);
        _acb_vec_add(new_z, new_z, z_mid, g, hprec);

        acb_theta_ql_all(all_val + k * n2, new_z, tau_mid, 0, hprec);
    }

    /* Make finite differences */
    for (k = 0; k < n2; k++)
    {
        for (j = 0; j < n_pow(b, g); j++)
        {
            acb_set(&val[j], &all_val[j * n2 + k]);
        }
        acb_theta_jet_ql_finite_diff(jet, eps, err, rho, val, ord, g, hprec);
        _acb_vec_set(dth + k * nb, jet, nb);
    }

    /* Add error */
    acb_theta_jet_naive_all(dth_low, z, tau, ord + 2, lp);
    for (k = 0; k < n2; k++)
    {
        acb_theta_jet_error_bounds(err_vec, z, tau, dth_low + k * nb_low, ord, lp);
        for (j = 0; j < nb; j++)
        {
            acb_add_error_arb(&dth[k * nb + j], &err_vec[j]);
        }
    }

    arb_clear(c);
    arb_clear(rho);
    arb_clear(t);
    arf_clear(eps);
    arf_clear(err);
    arf_clear(e);
    acb_mat_clear(tau_mid);
    _acb_vec_clear(z_mid, g);
    _acb_vec_clear(zetas, b);
    _acb_vec_clear(new_z, g);
    _acb_vec_clear(all_val, n2 * n_pow(b, g));
    _acb_vec_clear(val, n_pow(b, g));
    _acb_vec_clear(jet, nb);
    _acb_vec_clear(dth_low, n2 * nb_low);
    _arb_vec_clear(err_vec, nb);
}

void
acb_theta_jet_ql_all(acb_ptr dth, acb_srcptr z, const acb_mat_t tau, slong ord, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n2 = 1 << (2 * g);
    slong nb = acb_theta_jet_nb(ord, g);
    acb_ptr aux, new_z;
    arb_ptr v, a;
    acb_t c;
    arb_t u;
    slong k;

    aux = _acb_vec_init(nb);
    new_z = _acb_vec_init(g);
    v = _arb_vec_init(g);
    a = _arb_vec_init(g);
    acb_init(c);
    arb_init(u);

    acb_theta_naive_reduce(v, new_z, a, c, u, z, 1, tau, prec);
    acb_theta_jet_ql_all_red(dth, new_z, tau, ord, prec);

    _acb_vec_scalar_mul(dth, dth, n2 * nb, c, prec);
    _arb_vec_neg(a, a, g);
    _arb_vec_scalar_mul_2exp_si(a, a, g, 1);
    acb_theta_jet_exp_pi_i(aux, a, ord, g, prec);
    for (k = 0; k < n2; k++)
    {
        acb_theta_jet_mul(dth + k * nb, dth + k * nb, aux, ord, g, prec);
    }

    _acb_vec_clear(aux, nb);
    _acb_vec_clear(new_z, g);
    _arb_vec_clear(v, g);
    _arb_vec_clear(a, g);
    acb_clear(c);
    arb_clear(u);
}
