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
acb_theta_ctx_z_clear(acb_theta_ctx_z_t ctx)
{
    slong g = ctx->g;

    _acb_vec_clear(acb_theta_ctx_z(ctx), g);
    _acb_vec_clear(acb_theta_ctx_exp_z(ctx), g);
    acb_clear(acb_theta_ctx_c(ctx));
    _arb_vec_clear(acb_theta_ctx_r(ctx), g);
    arb_clear(acb_theta_ctx_uinv(ctx));
    arb_clear(acb_theta_ctx_u(ctx));

    if (g > 1)
    {
        _acb_vec_clear(acb_theta_ctx_exp_2z(ctx), g);
        _acb_vec_clear(acb_theta_ctx_exp_z_inv(ctx), g);
        _acb_vec_clear(acb_theta_ctx_exp_2z_inv(ctx), g);
        _arb_vec_clear(acb_theta_ctx_v(ctx), g);
    }
}
