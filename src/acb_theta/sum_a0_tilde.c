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

void acb_theta_sum_a0_tilde(acb_ptr th, const acb_theta_ctx_z_struct * vec, slong nb,
    const acb_theta_ctx_tau_t ctx_tau, arb_srcptr distances, slong prec)
{
    slong g = ctx_tau->g;
    slong n = 1 << g;
    slong new_prec;
    slong j;

    FLINT_ASSERT(nb >= 0);
    if (nb == 0)
    {
        return;
    }

    if (g == 1)
    {
        acb_ptr res;

        res = _acb_vec_init(4);
        new_prec = FLINT_MAX(prec + acb_theta_agm_addprec(&distances[0]),
            prec + acb_theta_agm_addprec(&distances[1]));

        for (j = 0; j < nb; j++)
        {
            acb_modular_theta_sum(&res[0], &res[1], &res[2], &res[3],
                (&vec[j])->exp_z, (&vec[j])->is_real,
                acb_mat_entry(ctx_tau->exp_tau, 0, 0), 1, new_prec);
            acb_set(&th[2 * j], &res[2]);
            acb_mul(&th[2 * j + 1], &res[1],
                acb_mat_entry(ctx_tau->exp_tau_div_4, 0, 0), new_prec);
            _acb_vec_scalar_mul_arb(th + 2 * j, th + 2 * j, 2,
                &(&vec[j])->uinv, new_prec);
        }

        _acb_vec_clear(res, 4);
    }
    else
    {
        /* Update the context for each a, call sum_00 with the right precision */
        acb_theta_ctx_z_struct * new_vec;
        acb_ptr res, cs;
        ulong a;

        new_vec = acb_theta_ctx_z_vec_init(nb, g);
        res = _acb_vec_init(nb);
        cs = _acb_vec_init(nb);

        for (a = 0; a < n; a++)
        {
            for (j = 0; j < nb; j++)
            {
                acb_theta_ctx_z_shift_a0(&new_vec[j], &cs[j], &vec[j], ctx_tau, a, prec);
            }
            new_prec = prec + acb_theta_agm_addprec(&distances[a]);
            acb_theta_sum_00(res, new_vec, nb, ctx_tau, new_prec);
            for (j = 0; j < nb; j++)
            {
                acb_mul(&res[j], &res[j], &cs[j], prec);
                acb_mul_arb(&th[n * j + a], &res[j], &(&vec[j])->uinv, prec);
            }
        }

        acb_theta_ctx_z_vec_clear(new_vec, nb);
        _acb_vec_clear(res, nb);
        _acb_vec_clear(cs, nb);
    }
}
