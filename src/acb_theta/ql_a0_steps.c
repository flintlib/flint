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

static void
acb_theta_ql_a0_step(acb_ptr r, acb_srcptr roots, arb_srcptr dist, arb_srcptr dist0,
    slong k, slong nb_steps, int has_t, int has_z, slong g, slong prec)
{
    slong n = 1 << g;
    arb_ptr d, d0;
    acb_ptr next;
    acb_ptr rts;
    slong nb_t = (has_t ? 3 : 1);
    slong nb_r = (has_t ? 2 : 1);
    slong nb_z = (has_z ? 2 : 1);
    slong j;

    d = _arb_vec_init(n);
    d0 = _arb_vec_init(n);
    next = _acb_vec_init(nb_z * nb_t * n);
    rts = _acb_vec_init(nb_r * nb_z * n);

    _arb_vec_scalar_mul_2exp_si(d, dist, n, k + 1);
    _arb_vec_scalar_mul_2exp_si(d0, dist0, n, k + 1);
    for (j = 0; j < nb_z * nb_r; j++)
    {
        _acb_vec_set(rts + j * n, roots + j * nb_steps * n + k * n, n);
    }

    if (has_t)
    {
        acb_theta_ql_step_3(next, r, r, rts, d0, d0, g, prec);
        if (has_z)
        {
            acb_theta_ql_step_3(next + nb_t * n, r + nb_t * n, r,
                rts + nb_r * n, d, d0, g, prec);
        }
    }
    else
    {
        acb_theta_ql_step_1(next, r, r, rts, d0, d0, g, prec);
        if (has_z)
        {
            acb_theta_ql_step_1(next + nb_t * n, r + nb_t * n, r,
                rts + nb_t * n, d, d0, g, prec);
        }
    }
    _acb_vec_set(r, next, nb_z * nb_t * n);

    _arb_vec_clear(d, n);
    _arb_vec_clear(d0, n);
    _acb_vec_clear(next, nb_z * nb_t * n);
    _acb_vec_clear(rts, nb_r * nb_z * n);
}

int
acb_theta_ql_a0_steps(acb_ptr r, acb_srcptr t, acb_srcptr z, arb_srcptr dist,
    arb_srcptr dist0, const acb_mat_t tau, slong guard, slong prec,
    acb_theta_ql_worker_t worker)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    int has_t = !_acb_vec_is_zero(t, g);
    int has_z = !_acb_vec_is_zero(z, g);
    slong nb_t = (has_t ? 3 : 1);
    slong nb_r = (has_t ? 2 : 1);
    slong nb_z = (has_z ? 2 : 1);
    acb_mat_t w;
    arb_mat_t Yinv;
    arb_mat_t cho;
    acb_ptr x, u, roots;
    arb_ptr y, new_dist;
    acb_t f, c;
    slong d, nb_steps;
    slong k;
    int res = 1;

    acb_mat_init(w, g, g);
    arb_mat_init(Yinv, g, g);
    arb_mat_init(cho, g, g);
    x = _acb_vec_init(g);
    u = _acb_vec_init(g);
    y = _arb_vec_init(g);
    new_dist = _arb_vec_init(n);
    acb_init(f);
    acb_init(c);

    /* Get f = i y Y^{-1} y */
    acb_mat_get_imag(Yinv, tau);
    arb_mat_inv(Yinv, Yinv, prec);
    _acb_vec_get_imag(y, z, g);
    arb_mat_bilinear_form(acb_imagref(f), Yinv, y, y, prec);

    /* Get nb_steps and dimension in ql_a0_split */
    acb_theta_eld_cho(cho, tau, ACB_THETA_LOW_PREC);
    d = acb_theta_ql_split(cho);
    nb_steps = acb_theta_ql_nb_steps(cho, d, prec);

    /* flint_printf("(ql_a0_steps) d =  %wd, has_z = %wd, has_t = %wd, cho:\n", d, has_z, has_t);
       arb_mat_printd(cho, 5); */
    /* flint_printf("(ql_a0_steps) Using d = %wd, nb_steps = %wd\n", d, nb_steps); */

    roots = _acb_vec_init(nb_z * nb_r * n * nb_steps);

    /* Get roots */
    res = acb_theta_ql_roots(roots, t, x, dist0, tau, nb_steps, guard, prec);
    if (res && has_z)
    {
        res = acb_theta_ql_roots(roots + nb_r * n * nb_steps, t, z, dist, tau,
            nb_steps, guard, prec);
    }

    if (res)
    {
        /* Call a0_split at 0 */
        acb_mat_scalar_mul_2exp_si(w, tau, nb_steps);
        _arb_vec_scalar_mul_2exp_si(new_dist, dist0, n, nb_steps);

        /* flint_printf("(ql_a0_steps) distances near cusp:\n");
        _arb_vec_printn(new_dist, n, 5, 0);
        flint_printf("\n"); */

        _acb_vec_scalar_mul_2exp_si(u, t, g, nb_steps);
        res = acb_theta_ql_a0_split(r, u, x, new_dist, w, d, guard, prec, worker);

        /* flint_printf("(ql_a0_steps) result of a0_split:\n");
        _acb_vec_printd(r, n * nb_t, 10);
        flint_printf("\n");*/
    }
    if (res && has_z)
    {
        /* Call a0_split at z and rescale */
        _acb_vec_scalar_mul_2exp_si(x, z, g, nb_steps);
        _arb_vec_scalar_mul_2exp_si(new_dist, dist, n, nb_steps);

        /*flint_printf("(ql_a0_steps) distances near cusp:\n");
        _arb_vec_printn(new_dist, n, 5, 0);
        flint_printf("\n");*/

        res = acb_theta_ql_a0_split(r + nb_t * n, u, x, new_dist, w, d,
            guard, prec, worker);
        acb_mul_2exp_si(c, f, nb_steps);
        acb_exp_pi_i(c, c, prec);
        _acb_vec_scalar_mul(r + nb_t * n, r + nb_t * n, n * nb_t, c, prec);

        /*flint_printf("(ql_a0_steps) result of a0_split and rescale:\n");
        _acb_vec_printd(r, n * nb_t, 10);
        flint_printf("\n");*/
    }

    if (res)
    {
        for (k = nb_steps - 1; k >= 0; k--)
        {
            acb_theta_ql_a0_step(r, roots, dist, dist0, k, nb_steps, has_t, has_z, g, prec);
            /*flint_printf("after step %wd\n", k);
            _acb_vec_printd(r, nb_z * nb_t * n, 5);
            flint_printf("\n");*/
        }
    }

    if (res && has_z)
    {
        acb_neg(c, f);
        acb_exp_pi_i(c, c, prec);
        _acb_vec_scalar_mul(r + nb_t * n, r + nb_t * n, n * nb_t, c, prec);
    }

    acb_mat_clear(w);
    arb_mat_clear(Yinv);
    arb_mat_clear(cho);
    _acb_vec_clear(x, g);
    _acb_vec_clear(u, g);
    _arb_vec_clear(y, g);
    _arb_vec_clear(new_dist, n);
    _acb_vec_clear(roots, nb_z * nb_r * n * nb_steps);
    acb_clear(f);
    acb_clear(c);
    return res;
}
