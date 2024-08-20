/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "acb.h"
#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"

void
acb_theta_ctx_init(acb_theta_ctx_t ctx, slong nb, slong g)
{
    slong n = 1 << g;

    FLINT_ASSERT(g >= 1);
    FLINT_ASSERT(nb >= 0);

    acb_mat_init(acb_theta_ctx_tau(ctx), g, g);
    arb_mat_init(acb_theta_ctx_y(ctx), g, g);
    arb_mat_init(acb_theta_ctx_yinv(ctx), g, g);
    acb_mat_init(acb_theta_ctx_exp_tau_div_4(ctx), g, g);
    acb_mat_init(acb_theta_ctx_exp_tau_div_2(ctx), g, g);
    acb_mat_init(acb_theta_ctx_exp_tau(ctx), g, g);
    acb_theta_ctx_exp_zs(ctx) = _acb_vec_init(nb * g);
    acb_theta_ctx_exp_zs_inv(ctx) = _acb_vec_init(nb * g);
    acb_theta_ctx_exp_2zs(ctx) = _acb_vec_init(nb * g);
    acb_theta_ctx_exp_2zs_inv(ctx) = _acb_vec_init(nb * g);
    acb_theta_ctx_cs(ctx) = _acb_vec_init(nb);
    acb_theta_ctx_us(ctx) = _arb_vec_init(nb);
    acb_theta_ctx_as(ctx) = _arb_vec_init(nb * g);
    acb_theta_ctx_nb(ctx) = nb;

    if (g >= 2)
    {
        arb_mat_init(acb_theta_ctx_cho(ctx), g, g);
        arb_mat_init(acb_theta_ctx_choinv(ctx), g, g);
        acb_mat_init(acb_theta_ctx_exp_tau_inv(ctx), g, g);
        acb_theta_ctx_vs(ctx) = _arb_vec_init(nb * g);
        acb_theta_ctx_d0(ctx) = _arb_vec_init(n);
        acb_theta_ctx_d(ctx) = _arb_vec_init(n);
    }

    ctx->t_is_zero = -1;
    ctx->z_is_zero = -1;
    ctx->z_is_real = -1;
}
