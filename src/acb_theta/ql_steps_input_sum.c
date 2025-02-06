/*
    Copyright (C) 2025 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"

void
acb_theta_ql_steps_input_sum(acb_ptr th_init, acb_srcptr zs, slong nb,
    acb_srcptr t, const acb_mat_t tau, arb_srcptr distances, slong nb_steps,
    const slong * easy_steps, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_theta_ctx_tau_t ctx_tau;
    acb_theta_ctx_z_struct * aux;
    acb_theta_ctx_z_t ctxt;
    acb_mat_t new_tau;
    acb_ptr new_z;
    arb_ptr d;
    slong j;

    acb_theta_ctx_tau_init(ctx_tau, 1, g);
    aux = acb_theta_ctx_z_vec_init(3, g);
    acb_mat_init(new_tau, g, g);
    new_z = _acb_vec_init(g);
    d = _arb_vec_init(n);
    if (easy_steps[0] < nb_steps)
    {
        acb_theta_ctx_z_init(ctxt, g);
    }

    acb_mat_scalar_mul_2exp_si(new_tau, tau, nb_steps);
    acb_theta_ctx_tau_set(ctx_tau, new_tau, prec);
    if (easy_steps[0] < nb_steps)
    {
        _acb_vec_scalar_mul_2exp_si(new_z, t, g, nb_steps);
        acb_theta_ctx_z_set(ctxt, new_z, ctx_tau, prec);
    }

    for (j = 0; j < nb; j++)
    {
        _acb_vec_scalar_mul_2exp_si(new_z, zs + j * g, g, nb_steps);
        acb_theta_ctx_z_set(&aux[0], new_z, ctx_tau, prec);
        _arb_vec_scalar_mul_2exp_si(d, distances + j * n, n, nb_steps);
        if (easy_steps[j] == nb_steps)
        {
            acb_theta_sum(th_init + 3 * n * j, aux, 1, ctx_tau, d, 1, 0, 1, prec);
        }
        else
        {
            acb_theta_ctx_z_add_real(&aux[1], &aux[0], ctxt, prec);
            acb_theta_ctx_z_add_real(&aux[2], &aux[1], ctxt, prec);
            if (j == 0)
            {
                acb_theta_sum(th_init + 3 * n * j, aux, 3, ctx_tau, d, 1, 0, 1, prec);
            }
            else
            {
                acb_theta_sum(th_init + 3 * n * j + n, aux + 1, 2, ctx_tau, d, 1, 0, 1, prec);
            }
        }
    }

    acb_theta_ctx_tau_clear(ctx_tau);
    acb_theta_ctx_z_vec_clear(aux, 3);
    acb_mat_clear(new_tau);
    _acb_vec_clear(new_z, g);
    _arb_vec_clear(d, n);
    if (easy_steps[0] < nb_steps)
    {
        acb_theta_ctx_z_clear(ctxt);
    }
}

