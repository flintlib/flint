/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

slong acb_theta_ql_reduce(acb_ptr new_z, acb_t c, arb_t u, ulong* a1, acb_srcptr z,
    const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_theta_eld_t E;
    arb_mat_t C, C1;
    acb_mat_t tau0, tau1, x;
    arb_ptr v, t, u;
    acb_t f;
    arf_t R2, eps;
    arb_t b;
    slong* n;
    slong s;

    arb_mat_init(C, g, g);
    v = _arb_vec_init(g);
    acb_init(f);
    arf_init(R2);
    arf_init(eps);
    arb_init(b);

    acb_theta_eld_C(C, tau, prec);
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
        arb_mat_init(C1, g - s, g - s);
        acb_mat_init(tau0, s, s);
        acb_mat_init(tau1, g - s, g - s);
        acb_mat_init(x, s, g - s);
        t = _acb_vec_init(g - s);
        u = _acb_vec_init(g - s);
        n = flint_malloc((g - s) * sizeof(slong));

        acb_theta_ql_blocks(tau0, x, tau1, tau, s);
        acb_theta_eld_cho(C1, t1, prec);
        _arb_vec_scalar_mul_2exp_si(v, v, g, 1);
        arf_mul_2exp_si(R2, R2, 2);
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
            acb_theta_eld_points(n, E);
            *a1 = acb_theta_char_get_a(n, g - s);

            /* Update new_z and c */
            acb_theta_char_get_acb(t, a, g - s);
            acb_mat_vector_mul_col(v, x, t, prec);
            _acb_vec_add(new_z, new_z, v, s, prec);

            acb_mat_vector_mul_col(v, tau1, t, prec);
            _acb_vec_scalar_mul_2exp_si(u, new_z + s, g - s, 1);
            _acb_vec_add(v, v, u, g - s, prec);
            acb_dot(f, NULL, 0, t, 1, v, 1, g - s, prec);
            acb_exp_pi_i(f, f, prec);
            acb_mul(c, c, f, prec);
        }

        acb_theta_eld_clear(E);
        arb_mat_clear(C1);
        acb_mat_clear(tau0);
        acb_mat_clear(tau1);
        acb_mat_clear(x);
        _acb_vec_clear(t, g - s);
        _acb_vec_init(u, g - s);
        flint_free(n);
    }

    arb_mat_clear(C);
    _arb_vec_clear(v, g);
    acb_clear(f);
    arf_clear(R2);
    arf_clear(eps);
    arb_clear(b);
    return s;
}
