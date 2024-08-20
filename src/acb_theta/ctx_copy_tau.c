/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"

void
acb_theta_ctx_copy_tau(acb_theta_ctx_t new_ctx, const acb_theta_ctx_t ctx)
{
    slong g = acb_theta_ctx_g(ctx);
    FLINT_ASSERT(acb_theta_ctx_g(new_ctx) == g);

    acb_mat_set(acb_theta_ctx_tau(new_ctx), acb_theta_ctx_tau(ctx));
    arb_mat_set(acb_theta_ctx_y(new_ctx), acb_theta_ctx_y(ctx));
    arb_mat_set(acb_theta_ctx_yinv(new_ctx), acb_theta_ctx_yinv(ctx));
    acb_mat_set(acb_theta_ctx_exp_tau_div_4(new_ctx), acb_theta_ctx_exp_tau_div_4(ctx));
    acb_mat_set(acb_theta_ctx_exp_tau_div_2(new_ctx), acb_theta_ctx_exp_tau_div_2(ctx));
    acb_mat_set(acb_theta_ctx_exp_tau(new_ctx), acb_theta_ctx_exp_tau(ctx));
    if (g > 1)
    {
        arb_mat_set(acb_theta_ctx_cho(new_ctx), acb_theta_ctx_cho(ctx));
        arb_mat_set(acb_theta_ctx_choinv(new_ctx), acb_theta_ctx_choinv(ctx));
        acb_mat_set(acb_theta_ctx_exp_tau_inv(new_ctx), acb_theta_ctx_exp_tau_inv(ctx));
    }
}
