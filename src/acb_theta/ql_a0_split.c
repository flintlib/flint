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
        arb_mat_indeterminate(Yinv);
    }

    _acb_vec_get_imag(res, z, g);
    arb_mat_vector_mul_col(res, Yinv, res, prec);

    arb_mat_clear(Yinv);
}

static void
acb_theta_ql_a0_eld_points(slong** pts, slong* nb_pts, arb_ptr v,
    slong* fullprec, arf_t eps, arb_srcptr d, ulong a, arb_srcptr nctr,
    const arb_mat_t C, const arb_mat_t C1, slong prec)
{
    slong g = arb_mat_nrows(C);
    slong s = g - arb_mat_nrows(C1);
    slong n = 1 << g;
    slong nba = 1 << (g - s);
    slong lp = ACB_THETA_LOW_PREC;
    arb_t max_d;
    arf_t R2;
    acb_theta_eld_t E;
    slong k;

    acb_theta_eld_init(E, g - s, g - s);
    arf_init(R2);
    arb_init(max_d);

    /* Get offset */
    acb_theta_char_get_arb(v, a, g - s);
    _arb_vec_add(v, v, nctr + s, g - s, prec);
    arb_mat_vector_mul_col(v, C1, v, prec);

    /* Get R2 */
    arb_zero(max_d);
    for (k = a; k < n; k += nba)
    {
        arb_max(max_d, max_d, &d[k], lp);
    }
    *fullprec = prec + acb_theta_dist_addprec(max_d);
    acb_theta_naive_radius(R2, eps, C, 0, *fullprec);

    /* List points in ellipsoid */
    acb_theta_eld_fill(E, C1, R2, v);
    *nb_pts = acb_theta_eld_nb_pts(E);
    *pts = flint_malloc(acb_theta_eld_nb_pts(E) * (g - s) * sizeof(slong));
    acb_theta_eld_points(*pts, E);

    acb_theta_eld_clear(E);
    arf_clear(R2);
    arb_init(max_d);
}

static int
acb_theta_ql_a0_split_term(acb_ptr th, slong* pt, ulong a, acb_srcptr t, acb_srcptr z,
    arb_srcptr v, arb_srcptr d, arb_srcptr new_d0, const acb_mat_t tau0,
    const acb_mat_t star, const acb_mat_t tau1, const arb_mat_t C1, slong guard,
    slong prec, slong fullprec, acb_theta_ql_worker_t worker)
{
    slong s = acb_mat_nrows(tau0);
    slong g = s + acb_mat_nrows(tau1);
    slong lp = ACB_THETA_LOW_PREC;
    slong n = 1 << g;
    slong nba = 1 << (g - s);
    slong nbth = 1 << s;
    slong nbt = (_acb_vec_is_zero(t, g) ? 1 : 3);
    slong new_prec;
    acb_ptr u, w, new_z, new_th;
    acb_t f, c;
    arb_ptr new_d;
    arb_t orth, x;
    slong j, k;
    int res;

    u = _acb_vec_init(g - s);
    w = _acb_vec_init(g - s);
    new_z = _acb_vec_init(s);
    new_th = _acb_vec_init(2 * nbth * nbt);
    new_d = _arb_vec_init(nbth);
    acb_init(f);
    acb_init(c);
    arb_init(orth);
    arb_init(x);

    /* Set u to pt + a1/2 */
    acb_theta_char_get_acb(u, a, g - s);
    for (j = 0; j < g - s; j++)
    {
        acb_add_si(&u[j], &u[j], pt[j], prec);
    }

    /* Get new_z and cofactor at 0 */
    acb_mat_vector_mul_col(new_z, star, u, prec);
    _acb_vec_add(new_z, new_z, z, s, prec);
    acb_dot(f, NULL, 0, u, 1, z + s, 1, g - s, prec);
    acb_mul_2exp_si(f, f, 1);
    acb_mat_vector_mul_col(w, tau1, u, prec);
    acb_dot(f, f, 0, w, 1, u, 1, g - s, prec);

    /* Get new distances and relative precision */
    acb_theta_dist_a0(new_d, new_z, tau0, lp);
    acb_theta_dist_pt(orth, v, C1, pt, lp);
    new_prec = prec;
    for (j = 0; j < nbth; j++)
    {
        arb_sub(x, &d[a + j * nba], orth, lp);
        arb_sub(x, x, &new_d[j], lp);
        new_prec = FLINT_MIN(new_prec, prec + acb_theta_dist_addprec(x));
    }
    new_prec = FLINT_MAX(new_prec, lp);

    /* Call worker */
    res = worker(new_th, t, new_z, new_d0, new_d, tau0, guard, new_prec);

    if (!_acb_vec_is_zero(new_z, s))
    {
        /* We are only interested in the values at z */
        _acb_vec_set(new_th, new_th + nbth * nbt, nbth * nbt);
    }

    /* Rescale to set th; cofactor depends on t */
    for (k = 0; k < nbt; k++)
    {
        acb_dot(c, NULL, 0, u, 1, t + s, 1, g - s, prec);
        acb_mul_si(c, c, 2 * k, prec);
        acb_add(c, c, f, prec);
        acb_exp_pi_i(c, c, prec);
        _acb_vec_scalar_mul(new_th + k * nbth, new_th + k * nbth,
            nbth, c, prec);
        for (j = 0; j < nbth; j++)
        {
            acb_add(&th[k * n + j * nba + a], &th[k * n + j * nba + a],
                &new_th[k * nbth + j], fullprec);
        }
    }

    _acb_vec_clear(u, g - s);
    _acb_vec_clear(w, g - s);
    _acb_vec_clear(new_z, s);
    _acb_vec_clear(new_th, 2 * nbth * nbt);
    _arb_vec_clear(new_d, nbth);
    acb_clear(f);
    acb_clear(c);
    arb_clear(orth);
    arb_clear(x);
    return res;
}

int
acb_theta_ql_a0_split(acb_ptr th, acb_srcptr t, acb_srcptr z, arb_srcptr d,
    const acb_mat_t tau, slong s, slong guard, slong prec, acb_theta_ql_worker_t worker)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong nba = 1 << (g - s);
    slong nbth = 1 << s;
    slong nbt = (_acb_vec_is_zero(t, g) ? 1 : 3);
    slong lp = ACB_THETA_LOW_PREC;
    arb_mat_t C, C1;
    acb_mat_t tau0, star, tau1;
    arb_ptr v, nctr, new_d0;
    arf_t eps;
    slong* pts;
    slong fullprec, nb_pts;
    slong a, j, k;
    int res = 1;

    if (s <= 0 || s >= g)
    {
        flint_printf("ql_a0_split: Error (must have 1 < s < g)\n");
        flint_abort();
    }

    arb_mat_init(C, g, g);
    arb_mat_init(C1, g - s, g - s);
    acb_mat_init(tau0, s, s);
    acb_mat_init(star, s, g - s);
    acb_mat_init(tau1, g - s, g - s);
    v = _arb_vec_init(g - s);
    nctr = _arb_vec_init(g);
    new_d0 = _arb_vec_init(nbth);
    arf_init(eps);

    acb_theta_ql_blocks(tau0, star, tau1, tau, s);
    acb_theta_eld_cho(C, tau, prec);
    for (j = 0; j < g - s; j++)
    {
        for (k = j; k < g - s; k++)
        {
            arb_set(arb_mat_entry(C1, j, k), arb_mat_entry(C, j, k));
        }
    }
    acb_theta_dist_a0(new_d0, z, tau0, lp);
    acb_theta_eld_ncenter(nctr, z, tau, prec);

    _acb_vec_zero(th, n * nbt);
    for (a = 0; a < nba; a++)
    {
        /* Get offset, fullprec, error and list of points in ellipsoid */
        acb_theta_ql_a0_eld_points(&pts, &nb_pts, v, &fullprec, eps,
            d, a, nctr, C, C1, prec);

        /* Sum terms at each point using worker */
        for (k = 0; (k < nb_pts) && res; k++)
        {
            res = acb_theta_ql_a0_split_term(th, pts + k * (g - s), a, t, z,
                v, d, new_d0, tau0, star, tau1, C1, guard, prec, fullprec, worker);
        }

        /* Add error */
        for (k = 0; k < nbth; k++)
        {
            for (j = 0; j < nbt; j++)
            {
                acb_add_error_arf(&th[j * n + k * nba + a], eps);
            }
        }

        flint_free(pts);
    }

    arb_mat_clear(C);
    arb_mat_clear(C1);
    acb_mat_clear(tau0);
    acb_mat_clear(star);
    acb_mat_clear(tau1);
    _arb_vec_clear(v, g - s);
    _arb_vec_clear(nctr, g);
    _arb_vec_clear(new_d0, nbth);
    arf_clear(eps);
    return res;
}
