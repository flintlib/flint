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
worker_dim1(acb_ptr th, acb_srcptr v1, acb_srcptr v2, const slong* precs, slong len,
    const acb_t cofactor, const slong* coords, slong ord, slong g, slong prec, slong fullprec)
{
    slong n = 1 << g;
    acb_t s0, s1, add, sub;
    ulong b;
    slong dot;

    acb_init(s0);
    acb_init(s1);
    acb_init(add);
    acb_init(sub);

    /* Compute alternate sums to adjust signs */
    acb_dot(s0, NULL, 0, v1, 2, v2, 2, (len + 1) / 2, prec);
    acb_dot(s1, NULL, 0, v1 + 1, 2, v2 + 1, 2, len / 2, prec);
    acb_add(add, s0, s1, prec);
    acb_sub(sub, s0, s1, prec);
    acb_mul(add, add, cofactor, prec);
    acb_mul(sub, sub, cofactor, prec);

    for (b = 0; b < n; b++)
    {
        dot = acb_theta_char_dot_slong(b, coords, g) % 2;
        if ((b >> (g - 1)) && dot)
        {
            acb_sub(&th[b], &th[b], sub, fullprec);
        }
        else if ((b >> (g - 1)))
        {
            acb_add(&th[b], &th[b], sub, fullprec);
        }
        else if (dot)
        {
            acb_sub(&th[b], &th[b], add, fullprec);
        }
        else
        {
            acb_add(&th[b], &th[b], add, fullprec);
        }
    }

    acb_clear(s0);
    acb_clear(s1);
    acb_clear(add);
    acb_clear(sub);
}

static void
acb_theta_naive_0b_gen(acb_ptr th, acb_srcptr z, slong nb_z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_theta_eld_t E;
    acb_theta_precomp_t D;
    acb_ptr c;
    arb_ptr u;
    acb_ptr new_z;
    slong ord = 0;
    slong nb = 1 << g;
    slong k;

    acb_theta_eld_init(E, g, g);
    acb_theta_precomp_init(D, nb_z, g);
    c = _acb_vec_init(nb_z);
    u = _arb_vec_init(nb_z);
    new_z = _acb_vec_init(nb_z * g);

    acb_theta_naive_ellipsoid(E, new_z, c, u, ord, z, nb_z, tau, prec);
    prec = acb_theta_naive_fullprec(E, prec);
    acb_theta_precomp_set(D, new_z, tau, E, prec);

    for (k = 0; k < nb_z; k++)
    {
        acb_theta_naive_worker_new(&th[k * nb], nb, &c[k], &u[k], E, D, k,
            ord, prec, worker_dim1);
    }

    acb_theta_eld_clear(E);
    acb_theta_precomp_clear(D);
    _acb_vec_clear(c, nb_z);
    _arb_vec_clear(u, nb_z);
    _acb_vec_clear(new_z, nb_z * g);
}

void
acb_theta_naive_0b(acb_ptr th, acb_srcptr z, slong nb_z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_ptr res;
    slong k;

    if (g == 1)
    {
        res = _acb_vec_init(4);
        for (k = 0; k < nb_z; k++)
        {
            acb_modular_theta(&res[0], &res[1], &res[2], &res[3], z + k * g,
                acb_mat_entry(tau, 0, 0), prec);
            acb_set(&th[2 * k], &res[2]);
            acb_set(&th[2 * k + 1], &res[3]);
        }
        _acb_vec_clear(res, 4);
    }
    else
    {
        acb_theta_naive_0b_gen(th, z, nb_z, tau, prec);
    }
}
