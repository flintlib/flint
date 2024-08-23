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

void acb_theta_sum_a0(acb_ptr th, const acb_theta_ctx_z_struct * vec, slong nb,
    const acb_theta_ctx_tau_t ctx_tau, arb_srcptr distances, slong prec)
{
    slong g = acb_theta_ctx_g(ctx_tau);
    slong n = 1 << g;
    acb_ptr res;
    acb_theta_ctx_z_struct * new_vec;
    slong new_prec;
    slong a, j;

    FLINT_ASSERT(nb >= 0);
    if (nb == 0)
    {
        return;
    }

    if (g == 1)
    {
        res = _acb_vec_init(4);
        for (j = 0; j < nb; j++)
        {
            /* acb_modular_theta_sum takes shifted precisions into account */
            acb_modular_theta_sum(&res[0], &res[1], &res[2], &res[3],
                acb_theta_ctx_exp_z(&vec[j]), acb_theta_ctx_is_real(&vec[j]),
                acb_mat_entry(acb_theta_ctx_exp_tau(ctx_tau), 0, 0), 1, prec);
            acb_mul(&th[2 * j], &res[2], acb_theta_ctx_c(&vec[j]), prec);
            acb_mul(&th[2 * j + 1], &res[1], acb_theta_ctx_c(&vec[j]), prec);
            acb_mul(&th[2 * j + 1], &th[2 * j + 1],
                acb_mat_entry(acb_theta_ctx_exp_tau_div_4(ctx_tau), 0, 0), prec);
        }
        _acb_vec_clear(res, 4);
    }
    else
    {
        /* Update the context for each a, call sum_00 with the right precision */
        new_vec = acb_theta_ctx_z_vec_init(nb, g);
        res = _acb_vec_init(nb);

        for (a = 0; a < n; a++)
        {
            for (j = 0; j < nb; j++)
            {
                acb_theta_ctx_z_shift_a0(&new_vec[j], &vec[j], ctx_tau, a, prec);
            }
            new_prec = prec + acb_theta_dist_addprec(&distances[a]);
            acb_theta_sum_00(res, new_vec, nb, ctx_tau, new_prec);
            for (j = 0; j < nb; j++)
            {
                acb_set(&th[n * j + a], &res[j]);
            }
        }

        acb_theta_ctx_z_vec_clear(new_vec, nb);
        _acb_vec_clear(res, nb);
    }
}
