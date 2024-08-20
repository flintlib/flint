/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_mat.h"
#include "acb_modular.h"
#include "acb_theta.h"

void
acb_theta_sum_all(acb_ptr th, const acb_theta_ctx_t ctx, slong start,
    slong nb, int z_is_real, slong prec)
{
    slong g = acb_theta_ctx_g(ctx);
    slong n = 1 << g;
    acb_ptr res;
    acb_theta_ctx_t new_ctx;
    slong new_prec;
    slong a, b, j, dot;

    if (g == 1)
    {
        res = _acb_vec_init(4);
        for (j = 0; j < nb; j++)
        {
            /* acb_modular_theta_sum takes shifted precisions into account */
            acb_modular_theta_sum(&res[0], &res[1], &res[2], &res[3],
                &acb_theta_ctx_exp_zs(ctx)[start + j], 0,
                acb_mat_entry(acb_theta_ctx_exp_tau(ctx), 0, 0), 1, prec);
            acb_set(&th[4 * j], &res[2]);
            acb_set(&th[4 * j + 1], &res[3]);
            acb_set(&th[4 * j + 2], &res[1]);
            acb_neg(&th[4 * j + 3], &res[0]);
            _acb_vec_scalar_mul(th + 4 * j + 2, th + 4 * j + 2, 2,
                acb_mat_entry(acb_theta_ctx_exp_tau_div_4(ctx), 0, 0), prec);
            _acb_vec_scalar_mul(th + 4 * j, th + 4 * j, 4,
                &acb_theta_ctx_cs(ctx)[start + j], prec);
        }
        _acb_vec_clear(res, 4);
    }
    else
    {
        /* Update the context for each a, call sum_00 with the right precision */
        acb_theta_ctx_init(new_ctx, nb, g);
        res = _acb_vec_init(n * nb);

        acb_theta_ctx_copy_tau(new_ctx, ctx);
        for (a = 0; a < n; a++)
        {
            acb_theta_ctx_shift_z(new_ctx, ctx, start, nb, a, prec);
            new_prec = prec + acb_theta_dist_addprec(
                (z_is_real ? &acb_theta_ctx_d0(ctx)[a] : &acb_theta_ctx_d(ctx)[a]));
            acb_theta_sum_0b(res, new_ctx, new_prec);
            for (j = 0; j < nb; j++)
            {
                _acb_vec_set(th + n * n * j + n * a, res + n * j, n);
            }
            for (b = 0; b < n; b++)
            {
                dot = acb_theta_char_dot(a, b, g);
                for (j = 0; j < nb; j++)
                {
                    acb_mul_i_pow_si(&th[n * n * j + n * a + b], &th[n * n * j + n * a + b], dot);
                }
            }
        }

        _acb_vec_clear(res, n * nb);
        acb_theta_ctx_clear(new_ctx);
    }
}
