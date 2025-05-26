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

static void
acb_theta_ctx_tau_set_shift(acb_theta_ctx_tau_t ctx, slong a, slong prec)
{
    slong g = ctx->g;
    acb_t x;
    slong j, k;

    acb_init(x);

    for (j = 0; j < g; j++)
    {
        /* Set exp_tau_a */
        acb_one(x);
        for (k = 0; k < g; k++)
        {
            if (!acb_theta_char_bit(a, k, g))
            {
                continue;
            }
            if (k < j)
            {
                acb_mul(x, x, acb_mat_entry(ctx->exp_tau_div_2, k, j), prec);
            }
            else if (k == j)
            {
                acb_mul(x, x, acb_mat_entry(ctx->exp_tau, k, k), prec);
            }
            else
            {
                acb_mul(x, x, acb_mat_entry(ctx->exp_tau_div_2, j, k), prec);
            }
        }
        acb_set(&ctx->exp_tau_a[a * g + j], x);

        /* Set exp_tau_a_inv; recompute products to avoid inversions */
        acb_one(x);
        for (k = 0; k < g; k++)
        {
            if (!acb_theta_char_bit(a, k, g))
            {
                continue;
            }
            if (k < j)
            {
                acb_mul(x, x, acb_mat_entry(ctx->exp_tau_div_2_inv, k, j), prec);
            }
            else if (k == j)
            {
                acb_mul(x, x, acb_mat_entry(ctx->exp_tau_inv, k, k), prec);
            }
            else
            {
                acb_mul(x, x, acb_mat_entry(ctx->exp_tau_div_2_inv, j, k), prec);
            }
        }
        acb_set(&ctx->exp_tau_a_inv[a * g + j], x);
    }
    /* Set exp_a_tau_a_div_4 */
    acb_one(x);
    for (j = 0; j < g; j++)
    {
        if (!acb_theta_char_bit(a, j, g))
        {
            continue;
        }
        for (k = j; k < g; k++)
        {
            if (!acb_theta_char_bit(a, k, g))
            {
                continue;
            }
            acb_mul(x, x, acb_mat_entry(ctx->exp_tau_div_4, j, k), prec);
        }
    }
    acb_set(&ctx->exp_a_tau_a_div_4[a], x);

    acb_clear(x);
}

void
acb_theta_ctx_tau_set(acb_theta_ctx_tau_t ctx, const acb_mat_t tau, slong prec)
{
    slong g = ctx->g;
    slong n = 1 << g;
    acb_t x;
    slong j, k, a;
    int b;

    FLINT_ASSERT(g == acb_mat_nrows(tau));
    acb_init(x);

    /* Set exp_tau and inverses */
    for (j = 0; j < g; j++)
    {
        for (k = j; k < g; k++)
        {
            acb_set_round(x, acb_mat_entry(tau, j, k), prec);
            acb_mul_2exp_si(x, x, (k == j ? -2 : -1));
            acb_exp_pi_i(acb_mat_entry(ctx->exp_tau_div_4, j, k), x, prec);

            acb_sqr(acb_mat_entry(ctx->exp_tau_div_2, j, k),
                acb_mat_entry(ctx->exp_tau_div_4, j, k), prec);
            acb_sqr(acb_mat_entry(ctx->exp_tau, j, k),
                acb_mat_entry(ctx->exp_tau_div_2, j, k), prec);

            /* Diagonal entries of exp_tau_inv are needed to set exp_tau_a_inv */
            b = acb_is_real(acb_mat_entry(tau, j, k));
            acb_theta_ctx_exp_inv(acb_mat_entry(ctx->exp_tau_div_4_inv, j, k),
                acb_mat_entry(ctx->exp_tau_div_4, j, k), x, b, prec);
            acb_theta_ctx_sqr_inv(acb_mat_entry(ctx->exp_tau_div_2_inv, j, k),
                acb_mat_entry(ctx->exp_tau_div_4_inv, j, k),
                acb_mat_entry(ctx->exp_tau_div_2, j, k), b, prec);
            acb_theta_ctx_sqr_inv(acb_mat_entry(ctx->exp_tau_inv, j, k),
                acb_mat_entry(ctx->exp_tau_div_2_inv, j, k),
                acb_mat_entry(ctx->exp_tau, j, k), b, prec);
        }
    }

    /* Set cho, yinv */
    acb_siegel_cho_yinv(&ctx->cho, &ctx->yinv, tau, prec);

    /* Set exponentials for shifts */
    if (ctx->allow_shift)
    {
        for (a = 0; a < n; a++)
        {
            acb_theta_ctx_tau_set_shift(ctx, a, prec);
        }
    }

    acb_clear(x);
}
