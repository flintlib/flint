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
#include "acb_theta.h"

void
acb_theta_ctx_z_init(acb_theta_ctx_z_t ctx, slong g)
{
    FLINT_ASSERT(g >= 1);

    ctx->g = g;
    acb_theta_ctx_z(ctx) = _acb_vec_init(g);
    acb_theta_ctx_exp_z(ctx) = _acb_vec_init(g);
    acb_init(acb_theta_ctx_c(ctx));
    acb_theta_ctx_r(ctx) = _arb_vec_init(g);

    if (g > 1)
    {
        acb_theta_ctx_exp_2z(ctx) = _acb_vec_init(g);
        acb_theta_ctx_exp_z_inv(ctx) = _acb_vec_init(g);
        acb_theta_ctx_exp_2z_inv(ctx) = _acb_vec_init(g);
        acb_theta_ctx_v(ctx) = _arb_vec_init(g);
        arb_init(acb_theta_ctx_u(ctx));
    }
}
