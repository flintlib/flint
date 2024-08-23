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
acb_theta_ctx_tau_clear(acb_theta_ctx_tau_t ctx)
{
    slong g = acb_theta_ctx_g(ctx);
    slong n = 1 << g;

    acb_mat_clear(acb_theta_ctx_tau(ctx));
    arb_mat_clear(acb_theta_ctx_y(ctx));
    arb_mat_clear(acb_theta_ctx_yinv(ctx));
    acb_mat_clear(acb_theta_ctx_exp_tau_div_4(ctx));
    acb_mat_clear(acb_theta_ctx_exp_tau_div_2(ctx));
    acb_mat_clear(acb_theta_ctx_exp_tau(ctx));

    if (g > 1)
    {
        arb_mat_clear(acb_theta_ctx_cho(ctx));
        arb_mat_clear(acb_theta_ctx_choinv(ctx));
        acb_mat_clear(acb_theta_ctx_exp_tau_div_4_inv(ctx));
        acb_mat_clear(acb_theta_ctx_exp_tau_div_2_inv(ctx));
        acb_mat_clear(acb_theta_ctx_exp_tau_inv(ctx));
        _acb_vec_clear(ctx->exp_tau_a_div_2, n * g);
        _acb_vec_clear(ctx->exp_tau_a, n * g);
        _acb_vec_clear(ctx->exp_tau_a_div_2_inv, n * g);
        _acb_vec_clear(ctx->exp_tau_a_inv, n * g);
        _acb_vec_clear(ctx->exp_a_tau_a_div_4, n);
    }
}
