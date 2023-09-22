/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

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
            ACB_THETA_QL_SPLIT);
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
    arb_mat_t cho;
    slong split, nb_steps;
    acb_mat_t tau_mid;
    acb_ptr t_mid, z_mid;
    arb_ptr err;
    slong k, j;
    int res;

    arb_mat_init(cho, g, g);
    acb_mat_init(tau_mid, g, g);
    t_mid = _acb_vec_init(g);
    z_mid = _acb_vec_init(g);
    err = _arb_vec_init(n * n);

    acb_theta_eld_cho(cho, tau, ACB_THETA_LOW_PREC);
    split = acb_theta_ql_split(cho);
    nb_steps = acb_theta_ql_nb_steps(cho, split, prec);

    /* Expect precision loss of (guard + g) * nb_steps, so call ql_a0_steps at midpoint */
    acb_mat_get_mid(tau_mid, tau);
    for (k = 0; k < g; k++)
    {
        for (j = 0; j <= k; j++)
        {
            acb_set(acb_mat_entry(tau_mid, j, k), acb_mat_entry(tau_mid, k, j));
        }
    }
    for (k = 0; k < g; k++)
    {
        acb_get_mid(&z_mid[k], &z[k]);
        acb_get_mid(&t_mid[k], &t[k]);
    }

    res = acb_theta_ql_a0_steps(r, t_mid, z_mid, dist0, dist, tau_mid, nb_steps,
        split, guard, prec + nb_steps * (guard/16 + g), &acb_theta_ql_a0);

    /* Add error */
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
        acb_theta_jet_error_bounds(err, z_mid, tau, 0, ACB_THETA_LOW_PREC);
        for (j = 0; j < n; j++)
        {
            acb_add_error_arb(&r[k * n + j], &err[n * j]);
        }
    }

    arb_mat_clear(cho);
    acb_mat_clear(tau_mid);
    _acb_vec_clear(t_mid, g);
    _acb_vec_clear(z_mid, g);
    _arb_vec_clear(err, n * n);
    return res;
}
