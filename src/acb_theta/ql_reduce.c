/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

slong acb_theta_ql_reduce(acb_ptr new_z, acb_t c, arb_t u, slong* n1, acb_srcptr z,
    const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_theta_eld_t E;
    arb_mat_t C, C1;
    acb_mat_t tau0, tau1, x;
    acb_ptr t, w;
    arb_ptr v;
    acb_t f;
    arf_t R2, eps;
    arb_t b;
    slong s, k;

    arb_mat_init(C, g, g);
    v = _arb_vec_init(g);
    acb_init(f);
    arf_init(R2);
    arf_init(eps);
    arb_init(b);

    acb_theta_eld_cho(C, tau, prec);
    acb_theta_naive_radius(R2, eps, C, 0, prec);
    acb_theta_naive_reduce(v, new_z, c, u, z, 1, tau, C, prec);
    arb_mul_arf(u, u, eps, prec);

    arb_set_arf(b, R2);
    arb_sqrt(b, b, prec);
    arb_mul_2exp_si(b, b, 1);

    for (s = g; s > 0; s--)
    {
        if (!arb_gt(arb_mat_entry(C, s - 1, s - 1), b))
        {
            break;
        }
    }

    if (s < g)
    {
        /* Construct ellipsoid */
        acb_theta_eld_init(E, g - s, g - s);
        arb_mat_window_init(C1, C, s, s, g, g);
        acb_mat_window_init(tau0, tau, 0, 0, s, s);
        acb_mat_window_init(tau1, tau, s, s, g, g);
        acb_mat_window_init(x, tau, 0, s, s, g);
        t = _acb_vec_init(g);
        w = _acb_vec_init(g);

        arb_mat_scalar_mul_2exp_si(C1, C1, -1);
        acb_theta_eld_fill(E, C1, R2, v + s);

        if (acb_theta_eld_nb_pts(E) == 0)
        {
            s = -1;
        }
        else if (acb_theta_eld_nb_pts(E) > 1)
        {
            flint_printf("(ql_reduce) Error: several points\n");
            flint_abort();
        }
        else
        {
            acb_theta_eld_points(n1, E);
            flint_printf("(ql_reduce) last coord: %wd\n", n1[g-s-1]);

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

            flint_printf("(ql_reduce) c, f:\n");
            acb_printd(c, 10);
            flint_printf("\n");
            acb_printd(f, 10);
            flint_printf("\n");
            
            acb_mul(c, c, f, prec);
        }

        acb_theta_eld_clear(E);
        arb_mat_window_clear(C1);
        acb_mat_window_clear(tau0);
        acb_mat_window_clear(tau1);
        acb_mat_window_clear(x);
        _acb_vec_clear(t, g);
        _acb_vec_clear(w, g);
    }

    arb_mat_clear(C);
    _arb_vec_clear(v, g);
    acb_clear(f);
    arf_clear(R2);
    arf_clear(eps);
    arb_clear(b);
    return s;
}
