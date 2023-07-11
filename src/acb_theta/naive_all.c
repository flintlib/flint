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
worker_dim0(acb_ptr th, const acb_t term, slong * coords, slong g,
            ulong ab, slong ord, slong prec, slong fullprec)
{
    acb_t x;
    slong sgn;
    ulong a = acb_theta_char_a(coords, g);
    ulong b;

    acb_init(x);

    for (b = 0; b < n_pow(2, g); b++)
    {
        sgn = acb_theta_char_dot_slong(b, coords, g) / 2;

        acb_set(x, term);
        if (sgn == 1)
            acb_mul_onei(x, x);
        else if (sgn == 2)
            acb_neg(x, x);
        else if (sgn == 3)
            acb_div_onei(x, x);

        ab = (a << g) + b;
        acb_add(&th[ab], &th[ab], x, fullprec);
    }

    acb_clear(x);
}

void
acb_theta_naive_all(acb_ptr th, acb_srcptr z, slong nb_z, const acb_mat_t tau,
                    slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_theta_eld_t E;
    acb_theta_precomp_t D;
    acb_ptr c;
    arf_t eps;
    int all = 1;
    slong ord = 0;
    ulong ab = 0;
    slong nb = 1 << (2 * g);
    acb_mat_t tau_adj;
    acb_ptr z_adj;
    slong k;

    acb_theta_eld_init(E, g, g);
    acb_theta_precomp_init(D, nb_z, g);
    c = _acb_vec_init(nb_z);
    acb_mat_init(tau_adj, g, g);
    z_adj = _acb_vec_init(g * nb_z);
    arf_init(eps);

    arf_one(eps);
    arf_mul_2exp_si(eps, eps, -prec);
    acb_theta_naive_ellipsoid(E, c, z_adj, ab, all, ord, z, nb_z, tau, eps, prec);
    prec = acb_theta_naive_fullprec(E, prec);

    acb_mat_scalar_mul_2exp_si(tau_adj, tau, -2);
    _acb_vec_scalar_mul_2exp_si(z_adj, z_adj, g * nb_z, -1);
    acb_theta_precomp_set(D, z_adj, tau_adj, E, prec);

    for (k = 0; k < nb_z; k++)
    {
        acb_theta_naive_worker(&th[k * nb], nb, &c[k], eps, E, D, k, ab,
                               ord, prec, worker_dim0);
    }

    acb_theta_eld_clear(E);
    acb_theta_precomp_clear(D);
    _acb_vec_clear(c, nb_z);
    acb_mat_clear(tau_adj);
    _acb_vec_clear(z_adj, g * nb_z);
    arf_clear(eps);
}
