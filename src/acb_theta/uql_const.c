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
agm_direct(acb_ptr th, acb_srcptr roots, const acb_mat_t tau,
    slong nb_steps, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_mat_t w;
    acb_ptr x;
    slong k;

    acb_mat_init(w, g, g);
    x = _acb_vec_init(g);

    acb_mat_scalar_mul_2exp_si(w, tau, nb_steps);
    acb_theta_naive_a0(th, x, 1, w, prec);
    
    for (k = nb_steps - 1; k >= 0; k--)
    {
        acb_theta_agm_sqr(th, th, g, prec);
        _acb_vec_scalar_mul_2exp_si(th, th, n, g);
        acb_theta_agm_sqrt(th, th, r + k * n, n, prec);
    }

    _acb_vec_clear(x, g);
    acb_mat_clear(w);
}

static void
agm_aux(acb_ptr th, acb_srcptr roots, acb_srcptr t, const acb_mat_t tau,
    slong nb_steps, slong prec)
{    
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_mat_t w;
    acb_ptr x;
    acb_ptr cur, next;
    slong k, a;
    
    acb_mat_init(w, g, g);
    x = _acb_vec_init(3 * g);
    cur = _acb_vec_init(3 * n);
    next = _acb_vec_init(3 * n);

    acb_mat_scalar_mul_2exp_si(w, tau, nb_steps);
    _acb_vec_scalar_mul_2exp_si(x + g, t, g, nb_steps);
    _acb_vec_scalar_mul_2exp_si(x + 2 * g, t, g, nb_steps + 1);    
    
    for (k = nb_steps - 1; k >= 0; k--)
    {
        /* Duplication formulas with square roots */
        acb_theta_agm_mul(next + n, cur, cur + n, g, hprec);
        acb_theta_agm_mul(next + 2 * n, cur, cur + 2 * n, g, hprec);
        _acb_vec_scalar_mul_2exp_si(next + n, next + n, 2 * n, g);
        acb_theta_agm_sqrt(next + n, next + n, r + k * 2 * n, 2 * n, hprec);
        
        /* Duplication formula with divisions */
        acb_theta_agm_sqr(next, cur + n, g, hprec);
        _acb_vec_scalar_mul_2exp_si(next, next, n, g);
        for (a = 0; a < n; a++)
        {
            acb_div(&next[a], &next[a], &next[2 * n + a], hprec);
        }
        _acb_vec_set(cur, next, 3 * n);
    }
    _acb_vec_set(th, cur, n);

    acb_mat_clear(w);
    _acb_vec_clear(x, 3 * g);
    _acb_vec_clear(cur, 3 * n);
    _acb_vec_clear(next, 3 * n);    
}

void
acb_theta_uql_const(acb_ptr th, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong nb_steps = acb_theta_ql_nb_steps(tau, prec);
    acb_ptr t;
    acb_ptr r;
    acb_ptr x;
    slong hprec;

    t = _acb_vec_init(g);
    r = _acb_vec_init(nb_steps * 2 * n);
    x = _acb_vec_init((nb_z + 1) * g);

    _acb_vec_set(x + g, z, nb_z * g);
        
    hprec = acb_theta_ql_roots(r, t, 1, tau, nb_steps, prec);
    if (hprec >= 0)
    {
        agm_direct(th, r, tau, nb_steps, hprec);
    }
    else
    {        
        hprec = acb_theta_uql_roots(r, t, t, 1, tau, nb_steps, prec);
        if (hprec >= 0)
        {
            agm_aux(th, r, t, tau, nb_steps, prec);
        }
        else
        {
            _acb_vec_indeterminate(th, n);
        }
    }

    _acb_vec_clear(t, g);
    _acb_vec_clear(r, nb_steps * 2 * n);
    _acb_vec_clear(x, (nb_z + 1) * g);
}
