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

static void
worker_dim1(acb_ptr dth, acb_srcptr v1, acb_srcptr v2, const slong* precs, slong len,
    const acb_t cofactor, const slong* coords, slong ord, slong g, slong prec, slong fullprec)
{
    slong n = 1 << g;
    acb_ptr v3, aux;
    acb_t x;
    slong a0, a1, b;
    slong* dots;
    slong i, ind;

    v3 = _acb_vec_init(len);
    aux = _acb_vec_init(2 * 3 * n);
    dots = flint_malloc(n * sizeof(slong));
    acb_init(x);

    /* Precompute a0, a1, dots and multiplications */
    a0 = acb_theta_char_get_a(coords, g);
    a1 = a0 ^ (1 << (g - 1));
    for (b = 0; b < n; b++)
    {
        dots[b] = acb_theta_char_dot_slong(b, coords, g);
    }
    for (i = 0; i < len; i++)
    {
        acb_mul(&v3[i], &v1[i], &v2[i], precs[i]);
    }

    /* Main loop */
    for (i = 0; i < len; i++)
    {
        for (b = 0; b < n; b++)
        {
            acb_mul_powi(x, &v3[i], (dots[b] + i * (b >> (g - 1))) % 4);
            ind = 3 * n * (i % 2) + 3 * b;
            if ((dots[b] + i * (b >> (g - 1))) % 2 == 0)
            {
                acb_add(&aux[ind], &aux[ind], x, prec);
            }
            else
            {
                acb_addmul_si(&aux[ind + 1], x, coords[0] + i, prec);
                acb_addmul_si(&aux[ind + 2], x, coords[1], prec);
            }
        }
    }

    /* Multiply vector by cofactor and i*pi, then add to dth */
    _acb_vec_scalar_mul(aux, aux, 2 * 3 * n, cofactor, prec);
    acb_const_pi(x, prec);
    acb_mul_onei(x, x);
    for (i = 0; i < 2 * n; i++)
    {
        _acb_vec_scalar_mul(&aux[3 * i + 1], &aux[3 * i + 1], 2, x, prec);
    }
    for (b = 0; b < n; b++)
    {
        _acb_vec_add(&dth[3 * (n * a0 + b)], &dth[3 * (n * a0 + b)],
            &aux[3 * b], 3, fullprec);
        _acb_vec_add(&dth[3 * (n * a1 + b)], &dth[3 * (n * a1 + b)],
            &aux[3 * (n + b)], 3, fullprec);
    }

    _acb_vec_clear(v3, len);
    _acb_vec_clear(aux, 2 * 3 * n);
    flint_free(dots);
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

    acb_theta_naive_worker_new(dth, 3 * n2, c, u, E, D, 0, ord, prec, worker_dim1);

    acb_theta_eld_clear(E);
    acb_theta_precomp_clear(D);
    _acb_vec_clear(z, g);
    acb_clear(c);
    arb_clear(u);
    acb_mat_clear(new_tau);
}
