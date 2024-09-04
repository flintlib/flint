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
acb_theta_ctx_tau_init(acb_theta_ctx_tau_t ctx, int allow_shift, slong g)
{
    slong n = 1 << g;
    FLINT_ASSERT(g >= 1);

    ctx->g = g;
    ctx->allow_shift = (g > 1) && allow_shift;
    arb_mat_init(&ctx->yinv, g, g);
    acb_mat_init(ctx->exp_tau_div_4, g, g);
    acb_mat_init(ctx->exp_tau_div_2, g, g);
    acb_mat_init(ctx->exp_tau, g, g);

    if (g > 1)
    {
        arb_mat_init(&ctx->cho, g, g);
        acb_mat_init(ctx->exp_tau_div_4_inv, g, g);
        acb_mat_init(ctx->exp_tau_div_2_inv, g, g);
        acb_mat_init(ctx->exp_tau_inv, g, g);

        if (allow_shift)
        {
            ctx->exp_tau_a_div_2 = _acb_vec_init(g * n);
            ctx->exp_tau_a = _acb_vec_init(g * n);
            ctx->exp_tau_a_div_2_inv = _acb_vec_init(g * n);
            ctx->exp_tau_a_inv = _acb_vec_init(g * n);
            ctx->exp_a_tau_a_div_4 = _acb_vec_init(n);
        }
    }
}
