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
worker_dim0(acb_ptr dth, slong nb, const acb_t term, slong* coords, slong g,
    slong ord, slong prec, slong fullprec)
{
    slong n = 1 << g;
    acb_t x;
    ulong a, b, ab;

    acb_init(x);

    a = acb_theta_char_get_a(coords, g);

    for (b = 0; b < n; b++)
    {
        ab = (a << g) + b;
        acb_mul_powi(x, term, acb_theta_char_dot_slong(b, coords, g) % 4);

        if (acb_theta_char_is_even(ab, 2))
        {
            acb_add(&dth[3 * ab], &dth[3 * ab], x, fullprec);
        }
        else
        {
            acb_addmul_si(&dth[3 * ab + 1], x, coords[0], fullprec);
            acb_addmul_si(&dth[3 * ab + 2], x, coords[1], fullprec);
        }
    }

    acb_clear(x);
}

void
acb_theta_g2_jet_naive_1(acb_ptr dth, const acb_mat_t tau, slong prec)
{
    slong g = 2;
    slong n2 = 1 << (2 * g);
    slong ord = 1;
    acb_theta_eld_t E;
    acb_theta_precomp_t D;
    acb_ptr z;
    acb_t c;
    arb_t u;
    acb_mat_t new_tau;

    acb_theta_eld_init(E, g, g);
    acb_theta_precomp_init(D, 1, g);
    z = _acb_vec_init(g);
    acb_init(c);
    arb_init(u);
    acb_mat_init(new_tau, g, g);

    acb_mat_scalar_mul_2exp_si(new_tau, tau, -2);

    acb_theta_jet_ellipsoid(E, u, z, new_tau, ord, prec);
    prec = acb_theta_naive_fullprec(E, prec);
    acb_theta_precomp_set(D, z, new_tau, E, prec);
    acb_one(c);

    acb_theta_naive_worker(dth, 3 * n2, c, u, E, D, 0, ord, prec, worker_dim0);

    acb_theta_eld_clear(E);
    acb_theta_precomp_clear(D);
    _acb_vec_clear(z, g);
    acb_clear(c);
    arb_clear(u);
    acb_mat_clear(new_tau);
}
