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
acb_theta_eld_ncenter(arb_ptr res, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arb_mat_t Yinv;
    int b;

    arb_mat_init(Yinv, g, g);

    acb_mat_get_imag(Yinv, tau);
    b = arb_mat_inv(Yinv, Yinv, prec);
    if (!b)
    {
        flint_printf("acb_theta_eld_center: Error (impossible inverse)\n");
        flint_printf("\n");
    }

    _acb_vec_get_imag(res, z, g);
    arb_mat_vector_mul_col(res, Yinv, res, prec);

    arb_mat_clear(Yinv);
}

static void
acb_theta_ql_blocks(acb_mat_t t0, acb_mat_t x, acb_mat_t t1,
    const acb_mat_t tau, slong d)
{
    slong g = acb_mat_nrows(tau);
    slong j, k;

    for (j = 0; j < d; j++)
    {
        for (k = 0; k < d; k++)
        {
            acb_set(acb_mat_entry(t0, j, k), acb_mat_entry(tau, j, k));
        }
    }

    for (j = 0; j < d; j++)
    {
        for (k = 0; k < (g - d); k++)
        {
            acb_set(acb_mat_entry(x, j, k), acb_mat_entry(tau, j, k + d));
        }
    }

    for (j = 0; j < (g - d); j++)
    {
        for (k = 0; k < (g - d); k++)
        {
            acb_set(acb_mat_entry(t1, j, k), acb_mat_entry(tau, j + d, k + d));
        }
    }
}

static void
acb_theta_ql_a0_eld_points(slong** pts, slong* nb_pts, arb_ptr offset,
    slong* fullprec, arf_t eps, arb_srcptr dist, ulong a, arb_srcptr nctr,
    const arb_mat_t cho, const arb_mat_t cho1, slong prec)
{
    slong g = arb_mat_nrows(cho);
    slong d = g - arb_mat_nrows(cho1);
    slong n = 1 << g;
    slong nb_a = 1 << (g - d);
    slong lp = ACB_THETA_LOW_PREC;
    arb_t max_dist;
    arf_t R2;
    acb_theta_eld_t E;
    slong k;

    acb_theta_eld_init(E, g - d, g - d);
    arf_init(R2);
    arb_init(max_dist);

    /* Get offset */
    acb_theta_char_get_arb(offset, a, g - d);
    _arb_vec_add(offset, offset, nctr + d, g - d, prec);
    arb_mat_vector_mul_col(offset, cho1, offset, prec);

    /* Get R2 */
    arb_zero(max_dist);
    for (k = a; k < n; k += nb_a)
    {
        arb_max(max_dist, max_dist, &dist[k], lp);
    }
    *fullprec = prec + acb_theta_dist_addprec(max_dist);
    acb_theta_naive_radius(R2, eps, cho, 0, *fullprec);

    /* List points in ellipsoid */
    acb_theta_eld_fill(E, cho1, R2, offset, prec);
    *nb_pts = acb_theta_eld_nb_pts(E);
    *pts = flint_malloc(acb_theta_eld_nb_pts(E) * (g - d) * sizeof(slong));
    acb_theta_eld_points(*pts, E);

    acb_theta_eld_clear(E);
    arf_clear(R2);
    arb_init(max_dist);
}

static int
acb_theta_ql_a0_split_term(acb_ptr r, slong* pt, ulong a, acb_srcptr t, acb_srcptr z,
    arb_srcptr offset, arb_srcptr dist, arb_srcptr new_dist0, const acb_mat_t tau0,
    const acb_mat_t star, const acb_mat_t tau1, const arb_mat_t cho1, slong guard,
    slong prec, slong fullprec, acb_theta_ql_worker_t worker)
{
    slong d = acb_mat_nrows(tau0);
    slong g = d + acb_mat_nrows(tau1);
    slong lp = ACB_THETA_LOW_PREC;
    slong n = 1 << g;
    slong nb_a = 1 << (g - d);
    slong nb_th = 1 << d;
    slong nb_t = (_acb_vec_is_zero(t, g) ? 1 : 3);
    slong new_prec;
    acb_ptr v, w, new_z, new_th;
    acb_t f, c;
    arb_ptr new_dist;
    arb_t orth_dist, x;
    slong j, k;
    int res;

    v = _acb_vec_init(g - d);
    w = _acb_vec_init(g - d);
    new_z = _acb_vec_init(d);
    new_th = _acb_vec_init(2 * nb_th * nb_t);
    new_dist = _arb_vec_init(nb_th);
    acb_init(f);
    acb_init(c);
    arb_init(orth_dist);
    arb_init(x);

    /* Set v to pt + a1/2 */
    acb_theta_char_get_acb(v, a, g - d);
    for (j = 0; j < g - d; j++)
    {
        acb_add_si(&v[j], &v[j], pt[j], prec);
    }

    /* Get new_z and cofactor at 0 */
    acb_mat_vector_mul_col(new_z, star, v, prec);
    _acb_vec_add(new_z, new_z, z, d, prec);
    acb_dot(f, NULL, 0, v, 1, z + d, 1, g - d, prec);
    acb_mul_2exp_si(f, f, 1);
    acb_mat_vector_mul_col(w, tau1, v, prec);
    acb_dot(f, f, 0, w, 1, v, 1, g - d, prec);

    /* Get new distances and relative precision */
    acb_theta_dist_a0(new_dist, new_z, tau0, lp);
    acb_theta_dist_pt(orth_dist, offset, cho1, pt, lp);
    new_prec = prec;
    for (j = 0; j < nb_th; j++)
    {
        arb_sub(x, &dist[a + j * nb_a], orth_dist, lp);
        arb_sub(x, x, &new_dist[j], lp);
        new_prec = FLINT_MIN(new_prec, prec + acb_theta_dist_addprec(x));
    }
    new_prec = FLINT_MAX(new_prec, lp);

    /* Call worker */
    /* flint_printf("nb_t = %wd, nb_th = %wd, d = %wd\n", nb_t, nb_th, d);
       flint_printf("new_prec = %wd, prec = %wd, fullprec = %wd, new_z:\n",
       new_prec, prec, fullprec);
       _acb_vec_printd(new_z, d, 5);
       flint_printf("\n");
       flint_printf("new_dist, orth_dist: ");
       _arb_vec_printn(new_dist, nb_th, 5, 0);
       flint_printf("\n");
       arb_printd(orth_dist, 5);
       flint_printf("\n"); */

    res = worker(new_th, t, new_z, new_dist0, new_dist, tau0, guard, new_prec);
    if (!_acb_vec_is_zero(new_z, d))
    {
        /* We are only interested in the values at z */
        _acb_vec_set(new_th, new_th + nb_th * nb_t, nb_th * nb_t);
    }

    /* flint_printf("output from worker:\n");
       _acb_vec_printd(new_th, nb_th * nb_t, 5);
       flint_printf("\n"); */

    /* Rescale to set r; cofactor depends on t */
    for (k = 0; k < nb_t; k++)
    {
        acb_dot(c, NULL, 0, v, 1, t + d, 1, g - d, prec);
        acb_mul_si(c, c, 2 * k, prec);
        acb_add(c, c, f, prec);
        acb_exp_pi_i(c, c, prec);

        /* flint_printf("cofactor for k = %wd: ", k);
           acb_printd(c, 10);
           flint_printf("\n");*/

        _acb_vec_scalar_mul(new_th + k * nb_th, new_th + k * nb_th,
            nb_th, c, prec);
        for (j = 0; j < nb_th; j++)
        {
            acb_add(&r[k * n + j * nb_a + a], &r[k * n + j * nb_a + a],
                &new_th[k * nb_th + j], fullprec);
        }
    }

    _acb_vec_clear(v, g - d);
    _acb_vec_clear(w, g - d);
    _acb_vec_clear(new_z, d);
    _acb_vec_clear(new_th, 2 * nb_th * nb_t);
    _arb_vec_clear(new_dist, nb_th);
    acb_clear(f);
    acb_clear(c);
    arb_clear(orth_dist);
    arb_clear(x);
    return res;
}

int
acb_theta_ql_a0_split(acb_ptr r, acb_srcptr t, acb_srcptr z, arb_srcptr dist,
    const acb_mat_t tau, slong d, slong guard, slong prec, acb_theta_ql_worker_t worker)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong nb_a = 1 << (g - d);
    slong nb_th = 1 << d;
    slong nb_t = (_acb_vec_is_zero(t, g) ? 1 : 3);
    slong lp = ACB_THETA_LOW_PREC;
    arb_mat_t cho, cho1;
    acb_mat_t tau0, star, tau1;
    arb_ptr offset, nctr, new_dist0;
    arf_t eps;
    slong* pts;
    slong fullprec, nb_pts;
    slong a, j, k;
    int res = 1;

    if (d <= 0 || d >= g)
    {
        flint_printf("ql_a0_split: Error (must have 1 < d < g - 1)\n");
        flint_abort();
    }

    arb_mat_init(cho, g, g);
    arb_mat_init(cho1, g - d, g - d);
    acb_mat_init(tau0, d, d);
    acb_mat_init(star, d, g - d);
    acb_mat_init(tau1, g - d, g - d);
    offset = _arb_vec_init(g - d);
    nctr = _arb_vec_init(g);
    new_dist0 = _arb_vec_init(nb_th);
    arf_init(eps);

    acb_theta_ql_blocks(tau0, star, tau1, tau, d);
    acb_theta_eld_cho(cho, tau, prec);
    acb_theta_eld_cho(cho1, tau1, prec);
    acb_theta_dist_a0(new_dist0, z, tau0, lp);
    acb_theta_eld_ncenter(nctr, z, tau, prec);

    _acb_vec_zero(r, n * nb_t);
    for (a = 0; a < nb_a; a++)
    {
        /* Get offset, fullprec, error and list of points in ellipsoid */
        acb_theta_ql_a0_eld_points(&pts, &nb_pts, offset, &fullprec, eps,
            dist, a, nctr, cho, cho1, prec);

        /* Sum terms at each point using worker */
        for (k = 0; (k < nb_pts) && res; k++)
        {
            res = acb_theta_ql_a0_split_term(r, pts + k * (g - d), a, t, z,
                offset, dist, new_dist0, tau0, star, tau1, cho1, guard,
                prec, fullprec,worker);
        }

        /* Add error */
        for (k = 0; k < nb_th; k++)
        {
            for (j = 0; j < nb_t; j++)
            {
                acb_add_error_arf(&r[j * n + k * nb_a + a], eps);
            }
        }

        flint_free(pts);
    }

    arb_mat_clear(cho);
    arb_mat_clear(cho1);
    acb_mat_clear(tau0);
    acb_mat_clear(star);
    acb_mat_clear(tau1);
    _arb_vec_clear(offset, g - d);
    _arb_vec_clear(nctr, g);
    _arb_vec_clear(new_dist0, nb_th);
    arf_clear(eps);
    return res;
}
