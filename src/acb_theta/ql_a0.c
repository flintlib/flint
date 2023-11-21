/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"

static slong
acb_theta_ql_split(const arb_mat_t cho)
{
    slong g = arb_mat_nrows(cho);
    arb_t cmp;
    slong k;

    arb_init(cmp);

    for (k = g - 1; k >= 1; k--)
    {
        arb_mul_2exp_si(cmp, arb_mat_entry(cho, k - 1, k - 1),
            FLINT_MAX(1, 6 + k - 2 * g));
        if (arb_lt(cmp, arb_mat_entry(cho, k, k)))
        {
            break;
        }
    }

    arb_clear(cmp);
    return k;
}

int
acb_theta_ql_a0(acb_ptr r, acb_srcptr t, acb_srcptr z, arb_srcptr dist0,
    arb_srcptr dist, const acb_mat_t tau, slong guard, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    int has_t = !_acb_vec_is_zero(t, g);
    int has_z = !_acb_vec_is_zero(z, g);
    slong nbt = (has_t ? 3 : 1);
    slong nbz = (has_z ? 2 : 1);
    slong nb_der = acb_theta_jet_nb(2, g);
    arb_mat_t cho;
    slong split, nb_steps, padding, lp;
    acb_mat_t tau_mid;
    acb_ptr t_mid, z_mid, dth;
    arb_t err;
    arf_t e;
    slong k, j, a;
    int res;

    arb_mat_init(cho, g, g);
    acb_mat_init(tau_mid, g, g);
    t_mid = _acb_vec_init(g);
    z_mid = _acb_vec_init(g);
    dth = _acb_vec_init(nb_der);
    arb_init(err);
    arf_init(e);

    acb_siegel_cho(cho, tau, ACB_THETA_LOW_PREC);
    split = acb_theta_ql_split(cho);
    nb_steps = acb_theta_ql_a0_nb_steps(cho, split, prec);
    if (has_t || has_z)
    {
        nb_steps += 1; /* cf p-ql_a0_steps */
    }
    padding = nb_steps * (guard + g);
    arf_one(e);
    arf_mul_2exp_si(e, e, -prec - padding);

    /* Expect precision loss of (guard + g) * nb_steps, so call ql_a0_steps at midpoint */
    acb_mat_get_mid(tau_mid, tau);
    for (k = 0; k < g; k++)
    {
        for (j = 0; j <= k; j++)
        {
            acb_add_error_arf(acb_mat_entry(tau_mid, k, j), e);
            acb_set(acb_mat_entry(tau_mid, j, k), acb_mat_entry(tau_mid, k, j));
        }
    }
    for (k = 0; k < g; k++)
    {
        acb_get_mid(&z_mid[k], &z[k]);
        acb_get_mid(&t_mid[k], &t[k]);
        if (has_z)
        {
            acb_add_error_arf(&z_mid[k], e);
        }
        if (has_t)
        {
            arb_add_error_arf(acb_realref(&t_mid[k]), e);
        }
    }

    res = acb_theta_ql_a0_steps(r, t_mid, z_mid, dist0, dist, tau_mid, nb_steps,
        split, guard, prec + padding, &acb_theta_ql_a0);

    /* Add error, using z_mid as temp */
    for (k = 0; (k < nbz * nbt) && res; k++)
    {
        _acb_vec_zero(z_mid, g);
        if (has_t)
        {
            _acb_vec_scalar_mul_ui(t_mid, t, g, k % 3, prec);
            _acb_vec_add(z_mid, z_mid, t_mid, g, prec);
        }
        if (has_z && (k >= nbt))
        {
            _acb_vec_add(z_mid, z_mid, z, g, prec);
        }
        for (a = 0; a < n; a++)
        {
            if (has_z && (k >= nbt))
            {
                lp = FLINT_MAX(ACB_THETA_LOW_PREC, acb_theta_dist_addprec(&dist[a]));
            }
            else
            {
                lp = FLINT_MAX(ACB_THETA_LOW_PREC, acb_theta_dist_addprec(&dist0[a]));
            }
            acb_theta_jet_naive_fixed_ab(dth, a << g, z_mid, tau, 2, lp);
            acb_theta_jet_error_bounds(err, z_mid, tau, dth, 0, lp);
            acb_add_error_arb(&r[k * n + a], err);
        }
    }

    arb_mat_clear(cho);
    acb_mat_clear(tau_mid);
    _acb_vec_clear(t_mid, g);
    _acb_vec_clear(z_mid, g);
    _acb_vec_clear(dth, nb_der);
    arb_clear(err);
    arf_clear(e);
    return res;
}
