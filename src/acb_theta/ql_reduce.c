/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"

slong
acb_theta_ql_reduce(acb_ptr new_z, acb_t c, arb_t u, slong * n1, acb_srcptr z,
    const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_theta_eld_t E;
    arb_mat_t C, W, C1;
    acb_mat_t tau0, tau1, x;
    acb_ptr t, w;
    arb_ptr v, a;
    acb_t f;
    arf_t R2, eps;
    slong s, k;
    int r = 1;

    arb_mat_init(C, g, g);
    v = _arb_vec_init(g);
    a = _arb_vec_init(g);
    acb_init(f);
    arf_init(R2);
    arf_init(eps);

    acb_siegel_cho(C, tau, prec);
    acb_theta_naive_radius(R2, eps, C, 0, prec);
    acb_theta_naive_reduce(v, new_z, a, c, u, z, 1, tau, prec);
    arb_mul_arf(u, u, eps, prec);

    for (s = g; (s >= 1) && r; )
    {
        s--;
        acb_theta_eld_init(E, g - s, g - s);
        arb_mat_window_init(W, C, s, s, g, g);
        arb_mat_init(C1, g - s, g - s);
        arb_mat_set(C1, W);

        arb_mat_scalar_mul_2exp_si(C1, C1, -1);
        r = acb_theta_eld_set(E, C1, R2, v + s);
        r = r && (acb_theta_eld_nb_pts(E) <= 1);
        if (r && (acb_theta_eld_nb_pts(E) == 0))
        {
            s = -2;
        }

        acb_theta_eld_clear(E);
        arb_mat_window_clear(W);
        arb_mat_clear(C1);
    }
    s++;

    if ((s >= 0) && (s < g))
    {
        /* We know E has exactly one point */
        acb_theta_eld_init(E, g - s, g - s);
        arb_mat_window_init(W, C, s, s, g, g);
        arb_mat_init(C1, g - s, g - s);
        acb_mat_window_init(tau0, tau, 0, 0, s, s);
        acb_mat_window_init(tau1, tau, s, s, g, g);
        acb_mat_window_init(x, tau, 0, s, s, g);
        t = _acb_vec_init(g);
        w = _acb_vec_init(g);

        arb_mat_set(C1, W);
        arb_mat_scalar_mul_2exp_si(C1, C1, -1);
        acb_theta_eld_set(E, C1, R2, v + s);
        acb_theta_eld_points(n1, E);

        /* Update new_z and c */
        for (k = 0; k < g - s; k++)
        {
            acb_set_si(&t[k], n1[k]);
        }
        _acb_vec_scalar_mul_2exp_si(t, t, g - s, -1);
        acb_mat_vector_mul_col(w, x, t, prec);
        _acb_vec_add(new_z, new_z, w, s, prec);

        acb_mat_vector_mul_col(w, tau1, t, prec);
        _acb_vec_scalar_mul_2exp_si(w, w, g - s, -1);
        _acb_vec_add(w, w, new_z + s, g - s, prec);
        _acb_vec_scalar_mul_2exp_si(w, w, g - s, 1);
        acb_dot(f, NULL, 0, t, 1, w, 1, g - s, prec);
        acb_exp_pi_i(f, f, prec);
        acb_mul(c, c, f, prec);

        acb_theta_eld_clear(E);
        arb_mat_window_clear(W);
        arb_mat_clear(C1);
        acb_mat_window_clear(tau0);
        acb_mat_window_clear(tau1);
        acb_mat_window_clear(x);
        _acb_vec_clear(t, g);
        _acb_vec_clear(w, g);
    }

    if (!arb_mat_is_finite(C)) /* early abort in ql_all */
    {
        acb_indeterminate(c);
        arb_pos_inf(u);
        s = -1;
    }

    arb_mat_clear(C);
    _arb_vec_clear(v, g);
    _arb_vec_clear(a, g);
    acb_clear(f);
    arf_clear(R2);
    arf_clear(eps);
    return s;
}
