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
agm_direct(acb_ptr th, acb_srcptr roots, acb_srcptr z, slong nb_z,
    const acb_mat_t tau, slong nb_steps, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_mat_t w;
    acb_ptr x;
    acb_ptr cur;
    slong k, j;

    acb_mat_init(w, g, g);
    x = _acb_vec_init((nb_z + 1) * g);
    cur = _acb_vec_init((nb_z + 1) * n);

    acb_mat_scalar_mul_2exp_si(w, tau, nb_steps);
    _acb_vec_scalar_mul_2exp_si(x + g, z, nb_z * g, nb_steps);
    acb_theta_naive_a0(cur, x, nb_z + 1, w, prec);

    for (k = nb_steps - 1; k >= 0; k--)
    {
        for (j = 0; j < nb_z; j++)
        {
            acb_theta_agm_mul(cur + (j + 1) * n, cur, cur + (j + 1) * n, g, prec);
        }
        acb_theta_agm_sqr(cur, cur, g, prec);
        acb_theta_agm_sqrt(cur, cur, r + k * (nb_z + 1) * n, prec);
    }
    _acb_vec_set(th, cur + n, nb_z * n);

    acb_mat_clear(w);
    _acb_vec_clear(x, (nb_z + 1) * g);
    _acb_vec_clear(cur, (nb_z + 1) * n);
}

static void
agm_aux(acb_ptr th, acb_srcptr roots, acb_srcptr t, acb_srcptr z, slong nb_z,
    const acb_mat_t tau, slong nb_steps, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_mat_t w;
    acb_ptr x;
    acb_ptr cur;
    acb_ptr next;
    slong k, j, a;

    acb_mat_init(w, g, g);
    x = _acb_vec_init(3 * (nb_z + 1) * g);
    cur = _acb_vec_init(3 * (nb_z + 1) * n);
    next = _acb_vec_init(3 * (nb_z + 1) * n);
        
    acb_mat_scalar_mul_2exp_si(w, tau, nb_steps);
    _acb_vec_set(x + g, t, g);
    _acb_vec_scalar_mul_2exp_si(x + 2 * g, t, g, 1);
    for (k = 0; k < nb_z; k++)
    {
        _acb_vec_set(x + (3 * k + 3) * g, z + k * g, prec);
        _acb_vec_add(x + (3 * k + 4) * g, x + g, z + k * g, prec);
        _acb_vec_add(x + (3 * k + 5) * g, x + 2 * g, z + k * g, prec);
    }
    _acb_vec_scalar_mul_2exp_si(x, x, 3 * (nb_z + 1) * g, nb_steps);
    acb_theta_naive_a0(cur, x, 3 * (nb_z + 1), w, prec);

    for (k = nb_steps - 1; k >= 0; k--)
    {
        /* Duplication using square roots */
        for (j = 0; j < nb_z + 1; j++)
        {
            acb_theta_agm_mul(next + (3 * j + 1) * n, cur, cur + (3 * j + 1) * n, prec);
            acb_theta_agm_mul(next + (3 * j + 2) * n, cur, cur + (3 * j + 2) * n, prec);
            _acb_vec_scalar_mul_2exp_si(next + (3 * j + 1) * n, next + (3 * j + 1) * n, 2 * j);
            acb_theta_agm_sqrt(next + (3 * j + 1) * n, next + (3 * j + 1) * n,
                r + k * 2 * j * n, 2 * n, prec);
        }

        /* Duplication using divisions */
        for (j = 0; j < nb_z + 1; j++)
        {
            acb_theta_agm_sqr(next + 3 * j * n, cur + (3 * j + 1) * n, g, prec);
            for (a = 0; a < n; a++)
            {
                acb_div(&next[3 * j * n + a], &next[3 * j * n + a],
                    &next[(3 * j + 2) * n + a], hprec);
            }                
        }
        _acb_vec_set(cur, next, 3 * (nb_z + 1) * n);
    }
    
    for (j = 0; j < nb_z; j++)
    {
        _acb_vec_set(th + j * n, cur + 3 * j * n);
    }
    
    acb_mat_clear(w);
    _acb_vec_clear(x, 3 * (nb_z + 1) * g);
    _acb_vec_clear(cur, 3 * (nb_z + 1) * n);
    _acb_vec_clear(next, 3 * (nb_z + 1) * n);
}

void
acb_theta_uql(acb_ptr th, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong nb_steps = acb_theta_ql_nb_steps(tau, prec);
    acb_ptr t;
    acb_ptr r;
    slong hprec;
    
    t = _acb_vec_init(g);
    r = _acb_vec_init(nb_steps * 2 * (nb_z + 1) * n);

    hprec = acb_theta_ql_roots(r, x, nb_z + 1, tau, nb_steps, prec);
    if (hprec >= 0)
    {
        agm_direct(th, roots, z, nb_z, tau, nb_steps, prec);
    }
    else
    {
        hprec = acb_theta_uql_roots(r, t, x, nb_z + 1, tau, nb_steps, prec);
        if (hprec >= 0)
        {
            agm_aux(th, r, t, z, nb_z, tau_nb_steps, prec);
        }
        else
        {
            _acb_vec_indeterminate(th, nb_z * n);
        }
    }
    
    _acb_vec_clear(t, g);
    _acb_vec_clear(r, nb_steps * 2 * nb_z * n);
}
