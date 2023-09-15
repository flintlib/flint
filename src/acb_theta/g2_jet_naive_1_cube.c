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
    ulong rem;
    slong p = n_rootrem(&rem, nb / (3 * n * n), 3);
    acb_ptr zetas = dth + nb;
    int is_even;
    slong dot;

    acb_t x;
    ulong a, b, ab;
    slong ind, k, a11, a12, a22, n1, n2;

    acb_init(x);

    a = acb_theta_char_get_a(coords, g);
    n1 = coords[0];
    n2 = coords[1];
    flint_printf("(worker) got nb = %wd, found p = %wd\n", nb, p);

    for (b = 0; b < n; b++)
    {
        ab = (a << g) + b;
        dot = acb_theta_char_dot_slong(b, coords, g);
        is_even = acb_theta_char_is_even(ab, 2);

        for (k = 0; k < n_pow(p, 3); k++)
        {
            a11 = k % p;
            a12 = (k / p) % p;
            a22 = (k / n_pow(p, 2)) % p;
            ind = (n1 * n1 * a11 + n2 * n2 * a22 + 2 * n1 * n2 * a12 + 2 * p * dot) % (8 * p);
            acb_mul(x, term, &zetas[ind], prec);
            ind = 3 * n * n * k + 3 * ab;

            if (n1 == 0 && n2 == 0)
            {
                flint_printf("n1 = 0, n2 = 0, a = %wd, k = %wd, ind = %wd, adding term:\n",
                    a, k, ind);
                acb_printd(term, 5);
                flint_printf("\n");
            }
            if (is_even)
            {
                acb_add(&dth[ind], &dth[ind], x, fullprec);
            }
            else
            {
                acb_addmul_si(&dth[ind + 1], x, n1, fullprec);
                acb_addmul_si(&dth[ind + 2], x, n2, fullprec);
            }
        }
    }

    acb_clear(x);
}

void
acb_theta_g2_jet_naive_1_cube(acb_ptr dth, const acb_mat_t tau, slong p, slong prec)
{
    slong g = 2;
    slong n2 = 1 << (2 * g);
    slong nb = 3 * n_pow(p, 3) * n2;
    slong ord = 1;
    acb_theta_eld_t E;
    acb_theta_precomp_t D;
    acb_ptr z, aux;
    acb_t c;
    arb_t u;
    acb_mat_t new_tau;

    acb_theta_eld_init(E, g, g);
    acb_theta_precomp_init(D, 1, g);
    z = _acb_vec_init(g);
    aux = _acb_vec_init(nb + 8 * p); /* cheat and add (8p)th roots at the end */
    acb_init(c);
    arb_init(u);
    acb_mat_init(new_tau, g, g);

    acb_mat_scalar_mul_2exp_si(new_tau, tau, -2);
    _acb_vec_unit_roots(aux + nb, 8 * p, 8 * p, prec);

    acb_theta_jet_ellipsoid(E, u, z, new_tau, ord, prec);
    prec = acb_theta_naive_fullprec(E, prec);
    acb_theta_precomp_set(D, z, new_tau, E, prec);
    acb_one(c);

    acb_theta_naive_worker(aux, nb, c, u, E, D, 0, ord, prec, worker_dim0);
    _acb_vec_set(dth, aux, nb);

    acb_theta_eld_clear(E);
    acb_theta_precomp_clear(D);
    _acb_vec_clear(z, g);
    acb_clear(c);
    arb_clear(u);
    acb_mat_clear(new_tau);
}
