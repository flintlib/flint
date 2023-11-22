/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"
#include "acb_theta.h"

#define ACB_THETA_G2_JET_NAIVE_1_THRESHOLD 100

static void
worker(acb_ptr dth, acb_srcptr v1, acb_srcptr v2, const slong * precs, slong len,
    const acb_t cofactor, const slong * coords, slong ord, slong g, slong prec, slong fullprec)
{
    slong n = 1 << g;
    acb_ptr v3, aux, sums_1, sums_2, diffs;
    acb_t x;
    slong a0, a1, b;
    slong * dots;
    slong i, ind0, ind1;

    v3 = _acb_vec_init(len);
    aux = _acb_vec_init(2 * 3 * n);
    sums_1 = _acb_vec_init(4);
    sums_2 = _acb_vec_init(4);
    diffs = _acb_vec_init(8);
    dots = flint_malloc(n * sizeof(slong));
    acb_init(x);

    /* Precompute a0, a1, dots and multiplications */
    a0 = acb_theta_char_get_a(coords, g);
    a1 = a0 ^ (1 << (g - 1));
    for (b = 0; b < n; b++)
    {
        dots[b] = acb_theta_char_dot_slong(b, coords, g);
    }

    if (len > ACB_THETA_G2_JET_NAIVE_1_THRESHOLD)
    {
        /* Store multiplications in v3 */
        for (i = 0; i < len; i++)
        {
            acb_mul(&v3[i], &v1[i], &v2[i], precs[i]);
        }
        /* Main loop */
        for (i = 0; i < len; i++)
        {
            for (b = 0; b < n; b++)
            {
                acb_mul_i_pow_si(x, &v3[i], (dots[b] + i * (b >> (g - 1))) % 4);
                ind0 = 3 * n * (i % 2) + 3 * b;
                if ((dots[b] + i * (b >> (g - 1))) % 2 == 0)
                {
                    acb_add(&aux[ind0], &aux[ind0], x, prec);
                }
                else
                {
                    acb_addmul_si(&aux[ind0 + 1], x, coords[0] + i, prec);
                    acb_addmul_si(&aux[ind0 + 2], x, coords[1], prec);
                }
            }
        }
    }
    else
    {
        /* Compute dot products and sum them appropriately */
        for (i = 0; i < len; i++)
        {
            acb_mul_si(&v3[i], &v2[i], coords[0] + i, precs[i]);
        }
        for (i = 0; i < 4; i++)
        {
            acb_dot(&sums_1[i], NULL, 0, v1 + i, 4, v2 + i, 4, (len + 3 - i) / 4, prec);
            acb_dot(&sums_2[i], NULL, 0, v1 + i, 4, v3 + i, 4, (len + 3 - i) / 4, prec);
        }
        acb_add(&diffs[0], &sums_1[0], &sums_1[2], prec);
        acb_add(&diffs[1], &sums_1[1], &sums_1[3], prec);
        acb_sub(&diffs[2], &sums_1[0], &sums_1[2], prec);
        acb_sub(&diffs[3], &sums_1[1], &sums_1[3], prec);
        acb_add(&diffs[4], &sums_2[0], &sums_2[2], prec);
        acb_add(&diffs[5], &sums_2[1], &sums_2[3], prec);
        acb_sub(&diffs[6], &sums_2[0], &sums_2[2], prec);
        acb_sub(&diffs[7], &sums_2[1], &sums_2[3], prec);

        /* Loop over b */
        for (b = 0; b < n; b++)
        {
            ind0 = 3 * b;
            ind1 = 3 * (n + b);
            if ((b >> (g - 1)) == 0)
            {
                if (dots[b] % 2 == 0)
                {
                    /* All even */
                    acb_mul_i_pow_si(x, &diffs[0], dots[b]);
                    acb_add(&aux[ind0], &aux[ind0], x, prec);
                    acb_mul_i_pow_si(x, &diffs[1], dots[b]);
                    acb_add(&aux[ind1], &aux[ind1], x, prec);
                }
                else
                {
                    /* All odd; use v3 for derivative wrt z1 */
                    acb_mul_i_pow_si(x, &diffs[4], dots[b]);
                    acb_add(&aux[ind0 + 1], &aux[ind0 + 1], x, prec);
                    acb_mul_i_pow_si(x, &diffs[0], dots[b]);
                    acb_add(&aux[ind0 + 2], &aux[ind0 + 2], x, prec);
                    acb_mul_i_pow_si(x, &diffs[5], dots[b]);
                    acb_add(&aux[ind1 + 1], &aux[ind1 + 1], x, prec);
                    acb_mul_i_pow_si(x, &diffs[1], dots[b]);
                    acb_add(&aux[ind1 + 2], &aux[ind1 + 2], x, prec);
                }
            }
            else
            {
                /* Alternating, with different signs for a0 and a1 */
                if (dots[b] % 2 == 0)
                {
                    /* a0 even, a1 odd */
                    acb_mul_i_pow_si(x, &diffs[2], dots[b]);
                    acb_add(&aux[ind0], &aux[ind0], x, prec);
                    acb_mul_i_pow_si(x, &diffs[7], dots[b] + 1);
                    acb_add(&aux[ind1 + 1], &aux[ind1 + 1], x, prec);
                    acb_mul_i_pow_si(x, &diffs[3], dots[b] + 1);
                    acb_add(&aux[ind1 + 2], &aux[ind1 + 2], x, prec);
                }
                else
                {
                    /* a0 odd, a1 even */
                    acb_mul_i_pow_si(x, &diffs[3], dots[b] + 1);
                    acb_add(&aux[ind1], &aux[ind1], x, prec);
                    acb_mul_i_pow_si(x, &diffs[6], dots[b]);
                    acb_add(&aux[ind0 + 1], &aux[ind0 + 1], x, prec);
                    acb_mul_i_pow_si(x, &diffs[2], dots[b]);
                    acb_add(&aux[ind0 + 2], &aux[ind0 + 2], x, prec);
                }
            }
        }
        /* Multiply d/dz2 entries by coords[1] */
        for (i = 0; i < 2 * n; i++)
        {
            acb_mul_si(&aux[3 * i + 2], &aux[3 * i + 2], coords[1], prec);
        }
    }

    /* Multiply vector by cofactor and add to dth */
    _acb_vec_scalar_mul(aux, aux, 2 * 3 * n, cofactor, prec);
    for (b = 0; b < n; b++)
    {
        _acb_vec_add(&dth[3 * (n * a0 + b)], &dth[3 * (n * a0 + b)],
            &aux[3 * b], 3, fullprec);
        _acb_vec_add(&dth[3 * (n * a1 + b)], &dth[3 * (n * a1 + b)],
            &aux[3 * (n + b)], 3, fullprec);
    }

    _acb_vec_clear(v3, len);
    _acb_vec_clear(aux, 2 * 3 * n);
    _acb_vec_clear(sums_1, 4);
    _acb_vec_clear(sums_2, 4);
    _acb_vec_clear(diffs, 8);
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
    acb_mat_t new_tau;
    arb_mat_t C;
    arf_t R2, eps;
    acb_ptr z;
    arb_ptr v, a;
    acb_t c;
    arb_t u;
    slong k;
    int b;

    acb_theta_eld_init(E, g, g);
    acb_mat_init(new_tau, g, g);
    arb_mat_init(C, g, g);
    arf_init(R2);
    arf_init(eps);
    z = _acb_vec_init(g);
    v = _arb_vec_init(g);
    a = _arb_vec_init(g);
    acb_init(c);
    arb_init(u);

    acb_mat_scalar_mul_2exp_si(new_tau, tau, -2);
    acb_siegel_cho(C, new_tau, prec);

    acb_theta_naive_reduce(v, z, a, c, u, z, 1, new_tau, prec);
    acb_theta_jet_naive_radius(R2, eps, C, v, ord, prec);
    b = acb_theta_eld_set(E, C, R2, v);

    if (b)
    {
        acb_theta_naive_worker(dth, 3 * n2, z, 1, new_tau, E, ord, prec, worker);

        arb_mul_arf(u, u, eps, prec);
        for (k = 0; k < 3 * n2; k++)
        {
            acb_add_error_arb(&dth[k], u);
        }

        _arb_vec_scalar_mul_2exp_si(a, a, 2, 1);
        _arb_vec_neg(a, a, 2);
        for (k = 0; k < n2; k++)
        {
            acb_addmul_arb(&dth[3 * k + 1], &dth[3 * k], &a[0], prec);
            acb_addmul_arb(&dth[3 * k + 2], &dth[3 * k], &a[1], prec);
        }

        acb_const_pi(c, prec);
        acb_mul_onei(c, c);
        for (k = 0; k < n2; k++)
        {
            acb_mul(&dth[3 * k + 1], &dth[3 * k + 1], c, prec);
            acb_mul(&dth[3 * k + 2], &dth[3 * k + 2], c, prec);
        }

    }
    else
    {
        for (k = 0; k < 3 * n2; k++)
        {
            acb_indeterminate(&dth[k]);
        }
    }

    acb_theta_eld_clear(E);
    acb_mat_clear(new_tau);
    arb_mat_clear(C);
    arf_clear(R2);
    arf_clear(eps);
    _acb_vec_clear(z, g);
    _arb_vec_clear(v, g);
    _arb_vec_clear(a, g);
    acb_clear(c);
    arb_clear(u);
}
