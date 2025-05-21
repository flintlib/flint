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

/* We make a choice between direct summation or ql_jet_fd */
/* See p-ql_jet_fd */

void
acb_theta_ql_jet(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong ord, int all, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong guard = ACB_THETA_LOW_PREC;
    slong * pattern;
    acb_theta_ctx_tau_t ctx_tau;
    acb_theta_ctx_z_struct * vec;
    slong j;
    int use_sum = 1;

    pattern = flint_malloc(g * sizeof(slong));

    acb_theta_ql_nb_steps(pattern, tau, 0, prec);
    if (pattern[0] >= 8 + ord
        || (g >= 2 && pattern[1] >= 6 + ord)
        || (g >= 3 && pattern[2] >= 7)
        || ord == 0)
    {
        use_sum = 0;
    }

    if (use_sum)
    {
        acb_theta_ctx_tau_init(ctx_tau, 0, g);
        vec = acb_theta_ctx_z_vec_init(nb, g);

        acb_theta_ctx_tau_set(ctx_tau, tau, prec + guard);
        for (j = 0; j < nb; j++)
        {
            acb_theta_ctx_z_set(&vec[j], zs + j * g, ctx_tau, prec + guard);
        }

        acb_theta_sum_jet(th, vec, nb, ctx_tau, ord, 1, all, prec);

        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_z_vec_clear(vec, nb);
    }
    else
    {
        acb_theta_ql_jet_fd(th, zs, nb, tau, ord, all, prec);
    }

    flint_free(pattern);
}
