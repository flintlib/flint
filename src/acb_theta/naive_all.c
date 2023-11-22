/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"
#include "acb_modular.h"
#include "acb_theta.h"

static void
acb_theta_naive_all_gen(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_ptr all_zs, ata, v;
    acb_t c;
    slong a, b, d, k;

    all_zs = _acb_vec_init(g * n * nb);
    ata = _acb_vec_init(n);
    v = _acb_vec_init(g);
    acb_init(c);

    for (a = 0; a < n; a++)
    {
        acb_theta_char_get_acb(v, a, g);
        acb_mat_vector_mul_col(v, tau, v, prec);
        for (k = 0; k < nb; k++)
        {
            _acb_vec_add(all_zs + k * g * n + a * g, zs + k * g, v, g, prec);
        }
        acb_theta_char_dot_acb(&ata[a], a, v, g, prec);
    }

    acb_theta_naive_0b(th, all_zs, n * nb, tau, prec);

    for (a = 0; a < n; a++)
    {
        /* Factors depending on z, not on b */
        for (k = 0; k < nb; k++)
        {
            acb_theta_char_dot_acb(c, a, zs + k * g, g, prec);
            acb_mul_2exp_si(c, c, 1);
            acb_add(c, c, &ata[a], prec);
            acb_exp_pi_i(c, c, prec);
            _acb_vec_scalar_mul(th + k * n * n + a * n,
                th + k * n * n + a * n, n, c, prec);
        }
        /* Factors depending on b, not on z */
        for (b = 0; b < n; b++)
        {
            d = acb_theta_char_dot(a, b, g);
            for (k = 0; k < nb; k++)
            {
                acb_mul_i_pow_si(&th[k * n * n + a * n + b],
                    &th[k * n * n + a * n + b], d);
            }
        }
    }

    _acb_vec_clear(all_zs, g * n * nb);
    _acb_vec_clear(ata, n);
    _acb_vec_clear(v, g);
    acb_clear(c);
}

static void
acb_theta_naive_all_g1(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec)
{
    acb_t q4, q, w;
    acb_ptr res;
    int w_is_unit;
    slong k;

    acb_init(q4);
    acb_init(q);
    acb_init(w);
    res = _acb_vec_init(4);

    acb_mul_2exp_si(q4, acb_mat_entry(tau, 0, 0), -2);
    acb_exp_pi_i(q4, q4, prec);
    acb_pow_ui(q, q4, 4, prec);

    for (k = 0; k < nb; k++)
    {
        acb_exp_pi_i(w, &zs[k], prec);
        w_is_unit = arb_is_zero(acb_imagref(&zs[k]));
        acb_modular_theta_sum(&res[0], &res[1], &res[2], &res[3],
            w, w_is_unit, q, 1, prec);
        acb_set(&th[4 * k], &res[2]);
        acb_set(&th[4 * k + 1], &res[3]);
        acb_mul(&th[4 * k + 2], &res[1], q4, prec);
        acb_mul(&th[4 * k + 3], &res[0], q4, prec);
        acb_neg(&th[4 * k + 3], &th[4 * k + 3]);
    }

    acb_clear(q4);
    acb_clear(q);
    acb_clear(w);
    _acb_vec_clear(res, 4);
}

void
acb_theta_naive_all(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);

    if (g == 1)
    {
        acb_theta_naive_all_g1(th, zs, nb, tau, prec);
    }
    else
    {
        acb_theta_naive_all_gen(th, zs, nb, tau, prec);
    }
}

