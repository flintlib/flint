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
    slong g = acb_theta_ctx_g(ctx);
    slong n = 1 << g;
    slong j, k;

    acb_mat_scalar_mul_2exp_si(acb_theta_ctx_tau(ctx), acb_theta_ctx_tau(ctx), 1);
    arb_mat_scalar_mul_2exp_si(acb_theta_ctx_y(ctx), acb_theta_ctx_y(ctx), 1);
    arb_mat_scalar_mul_2exp_si(acb_theta_ctx_yinv(ctx), acb_theta_ctx_yinv(ctx), -1);
    /* Swap matrices around */
    FLINT_SWAP(acb_mat_struct, *acb_theta_ctx_exp_tau_div_4(ctx),
        *acb_theta_ctx_exp_tau_div_2(ctx));
    FLINT_SWAP(acb_mat_struct, *acb_theta_ctx_exp_tau_div_2(ctx),
        *acb_theta_ctx_exp_tau(ctx));
    /* Update exp_tau */
    for (j = 0; j < g; j++)
    {
        for (k = j; k < g; k++)
        {
            acb_sqr(acb_mat_entry(acb_theta_ctx_exp_tau(ctx), j, k),
                acb_mat_entry(acb_theta_ctx_exp_tau_div_2(ctx), j, k), prec);
        }
    }

    if (g > 1)
    {
        arb_t sqrt2;
        acb_ptr temp;

        arb_init(sqrt2);
        arb_set_si(sqrt2, 2);
        arb_sqrt(sqrt2, sqrt2, prec);

        arb_mat_scalar_mul_arb(acb_theta_ctx_cho(ctx), acb_theta_ctx_cho(ctx), sqrt2, prec);
        arb_mat_scalar_div_arb(acb_theta_ctx_choinv(ctx), acb_theta_ctx_choinv(ctx), sqrt2, prec);
        FLINT_SWAP(acb_mat_struct, *acb_theta_ctx_exp_tau_div_4_inv(ctx),
            *acb_theta_ctx_exp_tau_div_2_inv(ctx));
        FLINT_SWAP(acb_mat_struct, *acb_theta_ctx_exp_tau_div_2_inv(ctx),
            *acb_theta_ctx_exp_tau_inv(ctx));
        for (j = 0; j < g; j++)
        {
            for (k = j; k < g; k++)
            {
                acb_theta_ctx_sqr_inv(acb_mat_entry(acb_theta_ctx_exp_tau_inv(ctx), j, k),
                    acb_mat_entry(acb_theta_ctx_exp_tau_div_2_inv(ctx), j, k),
                    acb_mat_entry(acb_theta_ctx_exp_tau(ctx), j, k),
                    acb_is_real(acb_mat_entry(acb_theta_ctx_tau(ctx), j, k)), prec);
            }
        }

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

        arb_clear(sqrt2);
    }
}
