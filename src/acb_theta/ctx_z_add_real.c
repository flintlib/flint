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
acb_theta_ctx_z_add_real(acb_theta_ctx_z_t res, const acb_theta_ctx_z_t ctx,
    const acb_theta_ctx_z_t ctx_real, slong prec)
{
    slong g = ctx->g;
    slong j;

    FLINT_ASSERT(ctx_real->g == g);
    FLINT_ASSERT(res->g == g);

    /* Copy things */
    arb_set(&res->u, &ctx->u);
    arb_set(&res->uinv, &ctx->uinv);
    res->is_real = ctx->is_real;
    if (g > 1)
    {
        _arb_vec_set(res->v, ctx->v, g);
    }

    /* Exponentials */
    for (j = 0; j < g; j++)
    {
        acb_mul(&res->exp_z[j], &ctx->exp_z[j], &ctx_real->exp_z[j], prec);
        if (g > 1)
        {
            acb_mul(&res->exp_z_inv[j], &ctx->exp_z_inv[j], &ctx_real->exp_z_inv[j], prec);
            acb_sqr(&res->exp_2z[j], &res->exp_z[j], prec);
            acb_theta_ctx_sqr_inv(&res->exp_2z_inv[j], &res->exp_z_inv[j],
                &res->exp_2z[j], ctx->is_real, prec);
        }
    }
}
