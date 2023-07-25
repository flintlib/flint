/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

static int
acb_theta_ql_roots_1(acb_ptr r, acb_srcptr z, arb_srcptr dist,
    const acb_t f, const acb_mat_t tau, slong nb_steps, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_mat_t w;
    acb_ptr x;
    acb_t c;
    arb_t d;
    slong hprec;
    slong k, a;
    int res = 1;

    acb_mat_init(w, g, g);
    x = _acb_vec_init(g);
    acb_init(c);
    arb_init(d);

    for (k = 0; (k < nb_steps) && res; k++)
    {
        acb_mat_scalar_mul_2exp_si(w, tau, k);
        _acb_vec_scalar_mul_2exp_si(x, z, g, k);
        acb_mul_2exp_si(c, f, k);
        acb_exp_pi_i(c, c, prec);

        for (a = 0; a < n; a++)
        {
            arb_mul_2exp_si(d, &dist[a], k);
            hprec = prec + acb_theta_dist_addprec(d);
            acb_theta_naive_ind(&r[k * n + a], a << g, x, 1, w, hprec);

            if (acb_contains_zero(&r[k * n + a]))
            {
                res = 0;
                break;
            }
        }

        _acb_vec_scalar_mul(r + k * n, r + k * n, n, c, prec);
    }

    acb_mat_clear(w);
    _acb_vec_clear(x, g);
    acb_clear(c);
    arb_clear(d);
    return res;
}

static int
acb_theta_ql_roots_3(acb_ptr r, acb_srcptr t, acb_srcptr z, arb_srcptr dist,
    const acb_mat_t tau, slong nb_steps, slong guard, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_ptr x;
    acb_t f;
    slong k;
    int res = 1;

    x = _acb_vec_init(g);
    acb_init(f);

    acb_theta_ql_log_rescale(f, z, tau, prec);

    if (_acb_vec_is_zero(t, g))
    {
        res = acb_theta_ql_roots_1(r, z, dist, f, tau, nb_steps, guard);
    }
    else
    {
        for (k = 1; (k < 3) && res; k++)
        {
            _acb_vec_scalar_mul_ui(x, t, g, k, prec);
            _acb_vec_add(x, x, z, g, prec);
            res = acb_theta_ql_roots_1(r + (k - 1) * nb_steps * n, x, dist,
                f, tau, nb_steps, guard);
        }
    }

    _acb_vec_clear(x, g);
    acb_clear(f);
    return res;
}

int
acb_theta_ql_roots(acb_ptr r, acb_srcptr t, acb_srcptr z, arb_srcptr dist0,
    arb_srcptr dist, const acb_mat_t tau, slong nb_steps, slong guard, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    int has_z = _acb_vec_is_zero(z, g);
    int has_t = _acb_vec_is_zero(t, g);
    slong nb_r = (has_t ? 2 : 1);
    acb_ptr x;
    int res;

    x = _acb_vec_init(g);

    res = acb_theta_ql_roots_3(r, t, x, dist0, tau, nb_steps, guard, prec);
    if (res && has_z)
    {
        res = acb_theta_ql_roots_3(r + nb_r * n * nb_steps, t, z, dist, tau,
            nb_steps, guard, prec);
    }

    _acb_vec_clear(x, g);
    return res;
}
