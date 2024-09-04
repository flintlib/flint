/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"
#include "acb_theta.h"

void acb_theta_00_notransform(acb_ptr th, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_theta_ctx_tau_t ctx_tau;
    acb_theta_ctx_z_struct * vec;
    slong j;

    acb_theta_ctx_tau_init(ctx_tau, 0, g);
    vec = acb_theta_ctx_z_vec_init(nb, g);

    acb_theta_ctx_tau_set(ctx_tau, tau, prec);
    for (j = 0; j < nb; j++)
    {
        acb_theta_ctx_z_set(&vec[j], zs + j * g, ctx_tau, prec);
    }

    acb_theta_sum_00(th, vec, nb, ctx_tau, prec);

    acb_theta_ctx_tau_clear(ctx_tau);
    acb_theta_ctx_z_vec_clear(vec, nb);
}
