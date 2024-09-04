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
acb_theta_ctx_tau_dupl(acb_theta_ctx_tau_t ctx, slong prec)
{
    slong g = ctx->g;
    slong n = 1 << g;
    slong j, k;

    arb_mat_scalar_mul_2exp_si(&ctx->yinv, &ctx->yinv, -1);

    /* Swap matrices around */
    FLINT_SWAP(acb_mat_struct, *ctx->exp_tau_div_4, *ctx->exp_tau_div_2);
    FLINT_SWAP(acb_mat_struct, *ctx->exp_tau_div_2, *ctx->exp_tau);

    /* Update exp_tau */
    for (j = 0; j < g; j++)
    {
        for (k = j; k < g; k++)
        {
            acb_sqr(acb_mat_entry(ctx->exp_tau, j, k),
                acb_mat_entry(ctx->exp_tau_div_2, j, k), prec);
        }
    }

    if (g > 1)
    {
        arb_t sqrt2;
        acb_ptr temp;

        arb_init(sqrt2);
        arb_set_si(sqrt2, 2);
        arb_sqrt(sqrt2, sqrt2, prec);

        arb_mat_scalar_mul_arb(&ctx->cho, &ctx->cho, sqrt2, prec);
        FLINT_SWAP(acb_mat_struct, *ctx->exp_tau_div_4_inv, *ctx->exp_tau_div_2_inv);
        FLINT_SWAP(acb_mat_struct, *ctx->exp_tau_div_2_inv, *ctx->exp_tau_inv);
        for (j = 0; j < g; j++)
        {
            for (k = j; k < g; k++)
            {
                acb_sqr(acb_mat_entry(ctx->exp_tau_inv, j, k),
                    acb_mat_entry(ctx->exp_tau_div_2_inv, j, k), prec);
            }
        }

        if (ctx->allow_shift)
        {
            temp = ctx->exp_tau_a_div_2;
            ctx->exp_tau_a_div_2 = ctx->exp_tau_a;
            ctx->exp_tau_a = temp;
            temp = ctx->exp_tau_a_div_2_inv;
            ctx->exp_tau_a_div_2_inv = ctx->exp_tau_a_inv;
            ctx->exp_tau_a_inv = temp;
            for (j = 0; j < n * g; j++)
            {
                acb_sqr(&ctx->exp_tau_a[j], &ctx->exp_tau_a_div_2[j], prec);
                acb_sqr(&ctx->exp_tau_a_inv[j], &ctx->exp_tau_a_div_2_inv[j], prec);
            }
            for (j = 0; j < n; j++)
            {
                acb_sqr(&ctx->exp_a_tau_a_div_4[j], &ctx->exp_a_tau_a_div_4[j], prec);
            }
        }

        arb_clear(sqrt2);
    }
}
