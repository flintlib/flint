/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "acb.h"
#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"

static void
acb_theta_ql_jet_all_sum(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong ord, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_theta_ctx_tau_t ctx_tau;
    acb_theta_ctx_z_struct * vec;
    slong j;

    acb_theta_ctx_tau_init(ctx_tau, g);
    vec = acb_theta_ctx_z_vec_init(nb, g);

    acb_theta_ctx_tau_set(ctx_tau, tau, prec);
    for (j = 0; j < nb; j++)
    {
        acb_theta_ctx_z_set(&vec[j], zs + j * g, ctx_tau, prec);
    }

    acb_theta_sum_jet_all(th, vec, nb, ctx_tau, ord, prec);

    acb_theta_ctx_tau_clear(ctx_tau);
    acb_theta_ctx_z_vec_clear(vec, nb);
}

static void
acb_theta_ql_jet_all_mid_err(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong ord, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n2 = 1 << (2 * g);
    slong b = ord + 1;
    slong hprec;
    slong lp = 8;
    slong nbth = acb_theta_jet_nb(ord, g);
    slong nb_low = acb_theta_jet_nb(ord + 2, g);
    slong nbpts = n_pow(b, g);
    arb_ptr c, rho, eps, err;
    arb_t t;
    acb_mat_t tau_mid;
    acb_ptr z_mid, zetas, new_zs;
    acb_ptr all_val, val, dth_low;
    acb_theta_ctx_tau_t ctx_tau;
    acb_theta_ctx_z_struct * vec;
    arb_ptr err_vec;
    slong k, kmod, j, l;

    c = _arb_vec_init(nb);
    rho = _arb_vec_init(nb);
    eps = _arb_vec_init(nb);
    err = _arb_vec_init(nb);
    arb_init(t);

    /* Get bounds and high precision, fail if too large */
    hprec = prec;
    for (l = 0; l < nb; l++)
    {
        acb_theta_jet_ql_bounds(&c[l], &rho[l], zs + l * g, tau, ord);
        acb_theta_jet_ql_radius(arb_midref(&eps[l]), arb_midref(&err[l]),
            &c[l], &rho[l], ord, g, prec);

        arb_log_base_ui(t, &eps[l], 2, lp);
        arb_neg(t, t);

        /* we expect that the second bound holds if finite, but still check */
        if (!arb_is_finite(t) || (arf_cmpabs_2exp_si(arb_midref(t), 20) > 0))
        {
            /* we could still compute for other zs, but this isn't supposed
               to happen anyway. */
            _acb_vec_indeterminate(th, nb * n2 * nbth);
            _arb_vec_clear(c, nb);
            _arb_vec_clear(rho, nb);
            _arb_vec_clear(eps, nb);
            _arb_vec_clear(err, nb);
            return;
        }

        hprec = FLINT_MAX(hprec, prec + ord * (arf_get_si(arb_midref(t), ARF_RND_CEIL) + g));
    }

    acb_mat_init(tau_mid, g, g);
    z_mid = _acb_vec_init(g);
    zetas = _acb_vec_init(b);
    new_zs = _acb_vec_init(g * nb * nbpts);
    all_val = _acb_vec_init(n2 * nb * nbpts);
    val = _acb_vec_init(nbpts);
    dth_low = _acb_vec_init(n2 * nb * nb_low);
    err_vec = _arb_vec_init(nbth);
    acb_theta_ctx_tau_init(ctx_tau, g);
    vec = acb_theta_ctx_z_vec_init(nb, g);

    /* Get midpoint of tau and zetas */
    for (j = 0; j < g; j++)
    {
        for (k = 0; k < g; k++)
        {
            acb_get_mid(acb_mat_entry(tau_mid, j, k), acb_mat_entry(tau, j, k));
            acb_set(acb_mat_entry(tau_mid, k, j), acb_mat_entry(tau_mid, j, k));
        }
    }
    _acb_vec_unit_roots(zetas, b, b, hprec);

    /* Fill in new_zs */
    for (l = 0; l < nb; l++)
    {
        for (j = 0; j < g; j++)
        {
            acb_get_mid(&z_mid[j], &zs[l * g + j]);
        }
        for (k = 0; k < nbpts; k++)
        {
            kmod = k;
            for (j = g - 1; j >= 0; j--)
            {
                acb_set(&new_zs[l * g * nbpts + k * g + j], &zetas[kmod % b]);
                kmod = kmod / b;
            }
            _acb_vec_scalar_mul_arb(new_zs + l * g * nbpts + k * g,
                new_zs + l * g * nbpts + k * g, g, &eps[l], hprec);
            _acb_vec_add(new_zs + l * g * nbpts + k * g,
                new_zs + l * g * nbpts + k * g, z_mid, g, hprec);
        }
    }

    /* Call theta */
    acb_theta_all_notransform(all_val, new_zs, nb * nbpts, tau_mid, 0, hprec);

    /* flint_printf("(jet_all_notransform) new_zs (%wd vectors):\n", nb * nbpts);
    _acb_vec_printd(new_zs, g * nb * nbpts, 5);
    flint_printf("(jet_all_notransform) got values:\n");
    _acb_vec_printd(all_val, nb * nbpts * n2, 5); */

    /* Make finite differences */
    for (l = 0; l < nb; l++)
    {
        for (k = 0; k < n2; k++)
        {
            for (j = 0; j < nbpts; j++)
            {
                acb_set(&val[j], &all_val[l * n2 * nbpts + j * n2 + k]);
            }
            acb_theta_jet_ql_finite_diff(th + l * n2 * nbth + k * nbth,
                arb_midref(&eps[l]), arb_midref(&err[l]), &rho[l], val, ord, g, hprec);
        }
    }

    /* Attempt to get finite error bounds */
    k = 0;
    do
    {
        k++;
        lp *= 2;
        acb_theta_ctx_tau_set(ctx_tau, tau, lp);
        for (l = 0; l < nb; l++)
        {
            acb_theta_ctx_z_set(&vec[l], zs + l * g, ctx_tau, lp);
        }
        acb_theta_sum_jet_all(dth_low, vec, nb, ctx_tau, ord + 2, lp);
    }
    while (!_acb_vec_is_finite(dth_low, nb * nb_low * n2) && k < 4);

    /* Add error */
    lp = ACB_THETA_LOW_PREC;
    for (l = 0; l < nb; l++)
    {
        for (k = 0; k < n2; k++)
        {
            acb_theta_jet_error_bounds(err_vec, zs + l * g, tau,
                dth_low + l * n2 * nb_low + k * nb_low, ord, lp);
            for (j = 0; j < nbth; j++)
            {
                acb_add_error_arb(&th[l * nbth * n2 + k * nbth + j], &err_vec[j]);
            }
        }
    }

    _arb_vec_clear(c, nb);
    _arb_vec_clear(rho, nb);
    _arb_vec_clear(eps, nb);
    _arb_vec_clear(err, nb);
    arb_clear(t);
    acb_mat_clear(tau_mid);
    _acb_vec_clear(z_mid, g);
    _acb_vec_clear(zetas, b);
    _acb_vec_clear(new_zs, nb * g * nbpts);
    _acb_vec_clear(all_val, nb * n2 * nbpts);
    _acb_vec_clear(val, nbpts);
    _acb_vec_clear(dth_low, nb * n2 * nb_low);
    _arb_vec_clear(err_vec, nbth);
    acb_theta_ctx_tau_clear(ctx_tau);
    acb_theta_ctx_z_vec_clear(vec, nb);
}

static int
acb_theta_jet_all_use_sum(const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong lp = ACB_THETA_LOW_PREC;
    arb_mat_t cho;
    slong * pattern;
    slong j;
    int b, res;

    arb_mat_init(cho, g, g);
    pattern = flint_malloc(g * sizeof(slong));

    acb_siegel_cho(cho, tau, lp);
    b = acb_theta_ql_nb_steps(pattern, cho, prec);

    /* Do not use sum only when at least 3 steps are needed */
    res = 1;
    if (b)
    {
        for (j = 0; j < g; j++)
        {
            if (pattern[j] > 2)
            {
                res = 0;
            }
        }
    }

    arb_mat_clear(cho);
    flint_free(pattern);
    return res;
}

void
acb_theta_jet_all_notransform(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong ord, slong prec)
{
    int use_sum = acb_theta_jet_all_use_sum(tau, prec);

    if (nb <= 0)
    {
        return;
    }

    if (use_sum)
    {
        acb_theta_ql_jet_all_sum(th, zs, nb, tau, ord, prec);
    }
    else
    {
        acb_theta_ql_jet_all_mid_err(th, zs, nb, tau, ord, prec);
    }
}
