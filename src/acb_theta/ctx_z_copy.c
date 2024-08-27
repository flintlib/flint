/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_theta.h"

void
acb_theta_ctx_z_copy(acb_theta_ctx_z_t res, const acb_theta_ctx_z_t ctx)
{
    slong g = ctx->g;

    FLINT_ASSERT(res->g == g);

    _acb_vec_set(acb_theta_ctx_z(res), acb_theta_ctx_z(ctx), g);
    _acb_vec_set(acb_theta_ctx_exp_z(res), acb_theta_ctx_exp_z(ctx), g);
    acb_set(acb_theta_ctx_c(res), acb_theta_ctx_c(ctx));
    _arb_vec_set(acb_theta_ctx_r(res), acb_theta_ctx_r(ctx), g);
    arb_set(acb_theta_ctx_u(res), acb_theta_ctx_u(ctx));
    arb_set(acb_theta_ctx_uinv(res), acb_theta_ctx_uinv(ctx));
    acb_theta_ctx_is_real(res) = acb_theta_ctx_is_real(ctx);

    if (g > 1)
    {
        _acb_vec_set(acb_theta_ctx_exp_2z(res), acb_theta_ctx_exp_2z(ctx), g);
        _acb_vec_set(acb_theta_ctx_exp_z_inv(res), acb_theta_ctx_exp_z_inv(ctx), g);
        _acb_vec_set(acb_theta_ctx_exp_2z_inv(res), acb_theta_ctx_exp_2z_inv(ctx), g);
        _arb_vec_set(acb_theta_ctx_v(res), acb_theta_ctx_v(ctx), g);
    }
}
