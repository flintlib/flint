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
acb_theta_sum_a0(acb_ptr th, const acb_theta_ctx_t ctx, slong start,
    slong nb, int z_is_real, slong prec)
{
    slong g = acb_theta_ctx_g(ctx);
    slong n = 1 << g;
    acb_ptr res;
    acb_theta_ctx_t new_ctx;
    slong new_prec;
    slong a, j;

    if (g == 1)
    {
        res = _acb_vec_init(4);
        for (j = 0; j < nb; j++)
        {
            /* acb_modular_theta_sum takes shifted precisions into account */
            acb_modular_theta_sum(&res[0], &res[1], &res[2], &res[3],
                &acb_theta_ctx_exp_zs(ctx)[start + j], 0,
                acb_mat_entry(acb_theta_ctx_exp_tau(ctx), 0, 0), 1, prec);
            acb_mul(&th[2 * j], &res[2], &acb_theta_ctx_cs(ctx)[start + j], prec);
            acb_mul(&th[2 * j + 1], &res[1], &acb_theta_ctx_cs(ctx)[start + j], prec);
            acb_mul(&th[2 * j + 1], &th[2 * j + 1],
                acb_mat_entry(acb_theta_ctx_exp_tau_div_4(ctx), 0, 0), prec);
        }
        _acb_vec_clear(res, 4);
    }
    else
    {
        /* Update the context for each a, call sum_00 with the right precision */
        acb_theta_ctx_init(new_ctx, nb, g);
        res = _acb_vec_init(nb);

        acb_theta_ctx_copy_tau(new_ctx, ctx);
        for (a = 0; a < n; a++)
        {
            acb_theta_ctx_shift_z(new_ctx, ctx, start, nb, a, prec);
            new_prec = prec + acb_theta_dist_addprec(
                (z_is_real ? &acb_theta_ctx_d0(ctx)[a] : &acb_theta_ctx_d(ctx)[a]));
            acb_theta_sum_00(res, new_ctx, new_prec);
            for (j = 0; j < nb; j++)
            {
                acb_set(&th[n * j + a], &res[j]);
            }
        }

        _acb_vec_clear(res, nb);
        acb_theta_ctx_clear(new_ctx);
    }
}
