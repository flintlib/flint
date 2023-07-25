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
acb_theta_ql_roots_one(acb_ptr r, acb_srcptr z, arb_srcptr dist,
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

int
acb_theta_ql_roots(acb_ptr r, acb_srcptr t, acb_srcptr z, arb_srcptr dist,
    const acb_mat_t tau, slong nb_steps, slong guard, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    arb_mat_t Yinv;
    arb_ptr y;
    acb_ptr x;
    acb_t f;
    slong k;
    int res = 1;

    arb_mat_init(Yinv, g, g);
    y = _arb_vec_init(g);
    x = _acb_vec_init(g);
    acb_init(f);

    /* Get i y Y^{-1} y */
    acb_mat_get_imag(Yinv, tau);
    arb_mat_inv(Yinv, Yinv, prec);
    _acb_vec_get_imag(y, z, g);
    arb_mat_bilinear_form(acb_imagref(f), Yinv, y, y, prec);

    if (_acb_vec_is_zero(t, g))
    {
        res = acb_theta_ql_roots_one(r, z, dist, f, tau, nb_steps, guard);
    }
    else
    {
        for (k = 1; (k < 3) && res; k++)
        {
            _acb_vec_scalar_mul_ui(x, t, g, k, prec);
            _acb_vec_add(x, x, z, g, prec);
            res = acb_theta_ql_roots_one(r + (k - 1) * nb_steps * n, x, dist,
                f, tau, nb_steps, guard);
        }
    }

    arb_mat_clear(Yinv);
    _arb_vec_clear(y, g);
    _acb_vec_clear(x, g);
    acb_clear(f);
    return res;
}
