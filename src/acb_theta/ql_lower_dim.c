/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"

static int
acb_theta_ql_eld_points(slong ** pts, slong * nb_pts, arb_ptr v,
    slong * fullprec, arf_t eps, arb_srcptr d, ulong a, arb_srcptr w,
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
    int res;

    acb_theta_eld_init(E, g - s, g - s);
    arf_init(R2);
    arb_init(max_d);

    /* Get offset */
    acb_theta_char_get_arb(v, a, g - s);
    _arb_vec_add(v, v, w + s, g - s, prec);
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
    res = acb_theta_eld_set(E, C1, R2, v);
    if (res)
    {
        *nb_pts = acb_theta_eld_nb_pts(E);
        *pts = flint_malloc(acb_theta_eld_nb_pts(E) * (g - s) * sizeof(slong));
        acb_theta_eld_points(*pts, E);
    }
    else
    {
        *nb_pts = 0;
        *pts = flint_malloc(0);
    }

    acb_theta_eld_clear(E);
    arf_clear(R2);
    arb_init(max_d);
    return res;
}

int
acb_theta_ql_lower_dim(acb_ptr * new_zs, acb_ptr * cofactors, slong ** pts,
    slong * nb, arf_t err, slong * fullprec, acb_srcptr z, const acb_mat_t tau,
    arb_srcptr distances, slong s, ulong a, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arb_mat_t C, C1, Yinv, Y0inv;
    acb_mat_t tau0, star, tau1;
    arb_ptr y, v, w, new_y, new_w;
    acb_ptr u, x;
    acb_t f;
    slong j, k;
    int res;

    FLINT_ASSERT(s >= 1 && s < g);

    arb_mat_init(C, g, g);
    arb_mat_init(Yinv, g, g);
    arb_mat_init(Y0inv, s, s);
    acb_mat_window_init(tau0, tau, 0, 0, s, s);
    acb_mat_window_init(star, tau, 0, s, s, g);
    acb_mat_window_init(tau1, tau, s, s, g, g);
    y = _arb_vec_init(g);
    v = _arb_vec_init(g - s);
    w = _arb_vec_init(g);
    u = _acb_vec_init(g - s);
    x = _acb_vec_init(g - s);
    new_y = _arb_vec_init(s);
    new_w = _arb_vec_init(s);
    acb_init(f);

    acb_siegel_yinv(Yinv, tau, prec);
    acb_siegel_yinv(Y0inv, tau0, prec);
    acb_siegel_cho(C, tau, prec);
    arb_mat_window_init(C1, C, s, s, g, g);
    _acb_vec_get_imag(y, z, g);
    arb_mat_vector_mul_col(w, Yinv, y, prec);

    res = acb_theta_ql_eld_points(pts, nb, v, fullprec,
        err, distances, a, w, C, C1, prec);
    *new_zs = _acb_vec_init((*nb) * s);
    *cofactors = _acb_vec_init(*nb);

    for (k = 0; k < *nb; k++)
    {
        /* Set u to pt + a1/2 */
        acb_theta_char_get_acb(u, a, g - s);
        for (j = 0; j < g - s; j++)
        {
            acb_add_si(&u[j], &u[j], (*pts)[k * (g - s) + j], prec);
        }

        /* Get new_z and cofactor at 0 */
        acb_mat_vector_mul_col(*new_zs + k * s, star, u, prec);
        _acb_vec_add(*new_zs + k * s, *new_zs + k * s, z, s, prec);
        acb_dot(f, NULL, 0, u, 1, z + s, 1, g - s, prec);
        acb_mul_2exp_si(f, f, 1);
        acb_mat_vector_mul_col(x, tau1, u, prec);
        acb_dot(f, f, 0, x, 1, u, 1, g - s, prec);

        arb_dot(acb_imagref(f), acb_imagref(f), 0, y, 1, w, 1, g, prec);
        _acb_vec_get_imag(new_y, *new_zs + k * s, s);
        arb_mat_vector_mul_col(new_w, Y0inv, new_y, prec);
        arb_dot(acb_imagref(f), acb_imagref(f), 1, new_y, 1, new_w, 1, s, prec);

        acb_exp_pi_i(*cofactors + k, f, prec);
    }

    arb_mat_clear(C);
    arb_mat_window_clear(C1);
    arb_mat_clear(Yinv);
    arb_mat_clear(Y0inv);
    acb_mat_window_clear(tau0);
    acb_mat_window_clear(star);
    acb_mat_window_clear(tau1);
    _arb_vec_clear(y, g);
    _arb_vec_clear(v, g - s);
    _arb_vec_clear(w, g);
    _acb_vec_clear(u, g - s);
    _acb_vec_clear(x, g - s);
    _arb_vec_clear(new_y, s);
    _arb_vec_clear(new_w, s);
    acb_clear(f);
    return res;
}
