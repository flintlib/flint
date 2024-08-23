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
acb_theta_ctx_tau_init(acb_theta_ctx_tau_t ctx, slong g)
{
    slong n = 1 << g;
    FLINT_ASSERT(g >= 1);

    acb_mat_init(acb_theta_ctx_tau(ctx), g, g);
    arb_mat_init(acb_theta_ctx_y(ctx), g, g);
    arb_mat_init(acb_theta_ctx_yinv(ctx), g, g);
    acb_mat_init(acb_theta_ctx_exp_tau_div_4(ctx), g, g);
    acb_mat_init(acb_theta_ctx_exp_tau_div_2(ctx), g, g);
    acb_mat_init(acb_theta_ctx_exp_tau(ctx), g, g);

    if (g > 1)
    {
        arb_mat_init(acb_theta_ctx_cho(ctx), g, g);
        arb_mat_init(acb_theta_ctx_choinv(ctx), g, g);
        acb_mat_init(acb_theta_ctx_exp_tau_div_4_inv(ctx), g, g);
        acb_mat_init(acb_theta_ctx_exp_tau_div_2_inv(ctx), g, g);
        acb_mat_init(acb_theta_ctx_exp_tau_inv(ctx), g, g);
        ctx->exp_tau_a_div_2 = _acb_vec_init(g * n);
        ctx->exp_tau_a = _acb_vec_init(g * n);
        ctx->exp_tau_a_div_2_inv = _acb_vec_init(g * n);
        ctx->exp_tau_a_inv = _acb_vec_init(g * n);
        ctx->exp_a_tau_a_div_4 = _acb_vec_init(n);
    }
}
