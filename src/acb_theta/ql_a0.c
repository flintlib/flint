/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

/* In the agm functions, guarantee nb_z >= 1 and first vector of z is zero */
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
    x = _acb_vec_init(nb_z * g);
    cur = _acb_vec_init(nb_z * n);

    acb_mat_scalar_mul_2exp_si(w, tau, nb_steps);
    _acb_vec_scalar_mul_2exp_si(x, z, nb_z * g, nb_steps);
    acb_theta_naive_a0(cur, x, nb_z, w, prec);

    for (k = nb_steps - 1; k >= 0; k--)
    {
        /* flint_printf("(ql_a0_direct) at step number %wd\n", k);
           _acb_vec_printd(cur, nb_z * n, 10);
           flint_printf("\n"); */
        
        for (j = 1; j < nb_z; j++)
        {
            acb_theta_agm_mul(cur + j * n, cur, cur + j * n, g, prec);
        }
        acb_theta_agm_sqr(cur, cur, g, prec);
        _acb_vec_scalar_mul_2exp_si(cur, cur, nb_z * n, g);
        
        /* flint_printf("(ql_a0) after duplication:\n");
           _acb_vec_printd(cur, nb_z * n, 10);
           flint_printf("\n");        
           flint_printf("(ql_a0) square roots:\n");
           _acb_vec_printd(roots + k * nb_z * n, nb_z * n, 10);
           flint_printf("\n"); */
        
        acb_theta_agm_sqrt(cur, cur, roots + k * nb_z * n, nb_z * n, prec);
    }
    _acb_vec_set(th, cur, nb_z * n);

    acb_mat_clear(w);
    _acb_vec_clear(x, nb_z * g);
    _acb_vec_clear(cur, nb_z * n);
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
    x = _acb_vec_init(3 * nb_z * g);
    cur = _acb_vec_init(3 * nb_z * n);
    next = _acb_vec_init(3 * nb_z * n);

    /* w = 2^k tau; x = 2^k (0, t, 2t, z1, z1 + t, z1 + 2t, ...) (since z0 = 0) */
    acb_mat_scalar_mul_2exp_si(w, tau, nb_steps);
    _acb_vec_set(x + g, t, g);
    _acb_vec_scalar_mul_2exp_si(x + 2 * g, t, g, 1);
    for (k = 1; k < nb_z; k++)
    {
        _acb_vec_set(x + (3 * k) * g, z + k * g, g);
        _acb_vec_add(x + (3 * k + 1) * g, x + g, z + k * g, g, prec);
        _acb_vec_add(x + (3 * k + 2) * g, x + 2 * g, z + k * g, g, prec);
    }
    _acb_vec_scalar_mul_2exp_si(x, x, 3 * nb_z * g, nb_steps);
    acb_theta_naive_a0(cur, x, 3 * nb_z, w, prec);

    for (k = nb_steps - 1; k >= 0; k--)
    {
        /*flint_printf("(ql_a0_aux) at step number %wd\n", k);
          _acb_vec_printd(cur, 3 * nb_z * n, 10);
          flint_printf("\n"); */
        
        /* Duplication using square roots for t, 2t, zi + t, zi + 2t */
        for (j = 0; j < nb_z; j++)
        {
            acb_theta_agm_mul(next + (3 * j + 1) * n, cur, cur + (3 * j + 1) * n, g, prec);
            acb_theta_agm_mul(next + (3 * j + 2) * n, cur, cur + (3 * j + 2) * n, g, prec);
            _acb_vec_scalar_mul_2exp_si(next + (3 * j + 1) * n,
                next + (3 * j + 1) * n, 2 * n, g);
        
            /* flint_printf("(ql_a0) after duplication:\n");
               _acb_vec_printd(next + (3 * j + 1) * n, 2 * n, 10);
               flint_printf("\n");        
               flint_printf("(ql_a0) square roots:\n");
               _acb_vec_printd(roots + k * 2 * nb_z * n + j * 2 * n, 2 * n, 10);
               flint_printf("\n"); */
            
            acb_theta_agm_sqrt(next + (3 * j + 1) * n, next + (3 * j + 1) * n,
                roots + k * 2 * nb_z * n + j * 2 * n, 2 * n, prec);
        }

        /* Duplication using divisions for 0 and zi */
        for (j = 0; j < nb_z; j++)
        {
            acb_theta_agm_mul(next + 3 * j * n, cur + (3 * j + 1) * n,
                cur + n, g, prec);
            _acb_vec_scalar_mul_2exp_si(next + 3 * j * n, next + 3 * j * n, n, g);
            for (a = 0; a < n; a++)
            {
                acb_div(&next[3 * j * n + a], &next[3 * j * n + a],
                    &next[(3 * j + 2) * n + a], prec);
            }                
        }
        _acb_vec_set(cur, next, 3 * nb_z * n);
        /*flint_printf("(ql_a0) after step number %wd\n", k);
        _acb_vec_printd(cur, 3 * nb_z * n, 10);
        flint_printf("\n");*/
    }
    
    for (j = 0; j < nb_z; j++)
    {
        _acb_vec_set(th + j * n, cur + 3 * j * n, n);
    }
    
    acb_mat_clear(w);
    _acb_vec_clear(x, 3 * nb_z * g);
    _acb_vec_clear(cur, 3 * nb_z * n);
    _acb_vec_clear(next, 3 * nb_z * n);
}

void
acb_theta_ql_a0(acb_ptr th, acb_srcptr z, slong nb_z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong nb_steps = acb_theta_ql_nb_steps(tau, prec);
    acb_ptr t;
    acb_ptr x;
    acb_ptr r;
    acb_ptr res;
    slong hprec;
    int has_zero = ((nb_z >= 1) && _acb_vec_is_zero(z, g));

    /* flint_printf("(acb_theta_ql_a0) g = %wd, prec = %wd, nb_steps = %wd, nb_z = %wd\n",
       g, prec, nb_steps, nb_z); */
    
    t = _acb_vec_init(g);
    x = _acb_vec_init((nb_z + 1) * g);
    r = _acb_vec_init(nb_steps * 2 * (nb_z + 1) * n);
    res = _acb_vec_init((nb_z + 1) * n);

    if (has_zero)
    {
        _acb_vec_set(x, z, nb_z * g);
        nb_z -= 1;
    }
    else
    {
        _acb_vec_set(x + g, z, nb_z * g);
    }
    
    hprec = acb_theta_ql_roots(r, x, nb_z + 1, tau, nb_steps, prec);
    if (hprec >= 0)
    {
        agm_direct(res, r, x, nb_z + 1, tau, nb_steps, hprec);
    }
    else
    {
        hprec = acb_theta_ql_roots_aux(r, t, x, nb_z + 1, tau, nb_steps, prec);
        if (hprec >= 0)
        {
            agm_aux(res, r, t, x, nb_z + 1, tau, nb_steps, hprec);
        }
        else
        {
            _acb_vec_indeterminate(res, (nb_z + 1) * n);
        }
    }

    if (has_zero)
    {
        _acb_vec_set(th, res, n * (nb_z + 1));
        nb_z += 1;
    }
    else
    {
        _acb_vec_set(th, res + n, n * nb_z);
    }
    
    _acb_vec_clear(t, g);
    _acb_vec_clear(x, (nb_z + 1) * g);
    _acb_vec_clear(r, nb_steps * 2 * (nb_z + 1) * n);
    _acb_vec_clear(res, (nb_z + 1) * n);
}
