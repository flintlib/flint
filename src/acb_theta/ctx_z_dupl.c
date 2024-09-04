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
acb_theta_ctx_z_dupl(acb_theta_ctx_z_t ctx, slong prec)
{
    slong g = ctx->g;
    acb_ptr temp;
    slong j;

    /* Compute exponentials, swapping vectors around if g > 1 */
    if (g == 1)
    {
        for (j = 0; j < g; j++)
        {
            acb_sqr(&ctx->exp_z[j], &ctx->exp_z[j], prec);
        }
    }
    else
    {
        temp = ctx->exp_z;
        ctx->exp_z = ctx->exp_2z;
        ctx->exp_2z = temp;
        temp = ctx->exp_z_inv;
        ctx->exp_z_inv = ctx->exp_2z_inv;
        ctx->exp_2z_inv = temp;
        for (j = 0; j < g; j++)
        {
            acb_sqr(&ctx->exp_2z[j], &ctx->exp_z[j], prec);
            acb_theta_ctx_sqr_inv(&ctx->exp_2z_inv[j], &ctx->exp_z_inv[j],
                &ctx->exp_2z[j], ctx->is_real, prec);
        }
    }

    /* Compute other quantities */
    arb_sqr(&ctx->u, &ctx->u, prec);
    arb_sqr(&ctx->uinv, &ctx->uinv, prec);

    if (g > 1)
    {
        arb_t sqrt2;

        arb_init(sqrt2);
        arb_set_si(sqrt2, 2);
        arb_sqrt(sqrt2, sqrt2, prec);

        _arb_vec_scalar_mul(ctx->v, ctx->v, g, sqrt2, prec);

        arb_clear(sqrt2);
    }
}
