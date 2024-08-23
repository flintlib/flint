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

static int
acb_theta_aj_is_zero(ulong a, slong j, slong g)
{
    return !((a >> (g - 1 - j)) & 1);
}

void
acb_theta_ctx_tau_set(acb_theta_ctx_tau_t ctx, const acb_mat_t tau, slong prec)
{
    slong g = acb_theta_ctx_g(ctx);
    slong n = 1 << g;
    acb_t x;
    slong j, k, a;
    int b;

    FLINT_ASSERT(g == acb_mat_nrows(tau));
    acb_init(x);

    /* Set tau, Y, exp_tau, and inverses */
    acb_mat_set(acb_theta_ctx_tau(ctx), tau);
    acb_mat_get_imag(acb_theta_ctx_y(ctx), tau);

    for (j = 0; j < g; j++)
    {
        for (k = j; k < g; k++)
        {
            acb_set_round(x, acb_mat_entry(tau, j, k), prec);
            acb_mul_2exp_si(x, x, (k == j ? -2 : -1));
            acb_exp_pi_i(acb_mat_entry(acb_theta_ctx_exp_tau_div_4(ctx), j, k), x, prec);

            acb_sqr(acb_mat_entry(acb_theta_ctx_exp_tau_div_2(ctx), j, k),
                acb_mat_entry(acb_theta_ctx_exp_tau_div_4(ctx), j, k), prec);
            acb_sqr(acb_mat_entry(acb_theta_ctx_exp_tau(ctx), j, k),
                acb_mat_entry(acb_theta_ctx_exp_tau_div_2(ctx), j, k), prec);

            if (g > 1)
            {
                /* Diagonal entries are also needed in shift_a0. */
                b = acb_is_real(acb_mat_entry(tau, j, k));
                acb_theta_ctx_exp_inv(acb_mat_entry(acb_theta_ctx_exp_tau_div_4_inv(ctx), j, k),
                    acb_mat_entry(acb_theta_ctx_exp_tau_div_4(ctx), j, k), x, b, prec);
                acb_theta_ctx_sqr_inv(acb_mat_entry(acb_theta_ctx_exp_tau_div_2_inv(ctx), j, k),
                    acb_mat_entry(acb_theta_ctx_exp_tau_div_4_inv(ctx), j, k),
                    acb_mat_entry(acb_theta_ctx_exp_tau_div_2(ctx), j, k), b, prec);
                acb_theta_ctx_sqr_inv(acb_mat_entry(acb_theta_ctx_exp_tau_inv(ctx), j, k),
                    acb_mat_entry(acb_theta_ctx_exp_tau_div_2_inv(ctx), j, k),
                    acb_mat_entry(acb_theta_ctx_exp_tau(ctx), j, k), b, prec);
            }
        }
    }

    /* Set C, Cinv, Yinv */
    if (g == 1)
    {
        arb_inv(arb_mat_entry(acb_theta_ctx_yinv(ctx), 0, 0),
            arb_mat_entry(acb_theta_ctx_y(ctx), 0, 0), prec);
    }
    else
    {
        arb_t pi;
        arb_init(pi);

        arb_const_pi(pi, prec);
        acb_siegel_cho(acb_theta_ctx_cho(ctx), tau, prec); /* has a factor sqrt(pi) */
        b = arb_mat_inv(acb_theta_ctx_choinv(ctx), acb_theta_ctx_cho(ctx), prec);
        if (!b)
        {
            arb_mat_indeterminate(acb_theta_ctx_choinv(ctx));
        }
        arb_mat_transpose(acb_theta_ctx_yinv(ctx), acb_theta_ctx_choinv(ctx));
        arb_mat_mul(acb_theta_ctx_yinv(ctx), acb_theta_ctx_choinv(ctx), acb_theta_ctx_yinv(ctx), prec);
        arb_mat_scalar_mul_arb(acb_theta_ctx_yinv(ctx), acb_theta_ctx_yinv(ctx), pi, prec);

        arb_clear(pi);
    }

    /* Set exponentials for shifts */
    if (g > 1)
    {
        for (a = 0; a < n; a++)
        {
            for (j = 0; j < g; j++)
            {
                acb_one(x);
                for (k = 0; k < g; k++)
                {
                    if (acb_theta_aj_is_zero(a, k, g))
                    {
                        continue;
                    }
                    if (k < j)
                    {
                        acb_mul(x, x, acb_mat_entry(acb_theta_ctx_exp_tau_div_4(ctx), k, j), prec);
                    }
                    else if (k == j)
                    {
                        acb_mul(x, x, acb_mat_entry(acb_theta_ctx_exp_tau_div_2(ctx), k, k), prec);
                    }
                    else
                    {
                        acb_mul(x, x, acb_mat_entry(acb_theta_ctx_exp_tau_div_4(ctx), j, k), prec);
                    }
                }
                acb_set(&acb_theta_ctx_exp_tau_a_div_2(ctx, a)[j], x);
                acb_sqr(&acb_theta_ctx_exp_tau_a(ctx, a)[j], x, prec);
                /* Recompute multiplications to avoid inversions. */
                acb_one(x);
                for (k = 0; k < g; k++)
                {
                    if (acb_theta_aj_is_zero(a, k, g))
                    {
                        continue;
                    }
                    if (k < j)
                    {
                        acb_mul(x, x, acb_mat_entry(acb_theta_ctx_exp_tau_div_4_inv(ctx), k, j), prec);
                    }
                    else if (k == j)
                    {
                        acb_mul(x, x, acb_mat_entry(acb_theta_ctx_exp_tau_div_2_inv(ctx), k, k), prec);
                    }
                    else
                    {
                        acb_mul(x, x, acb_mat_entry(acb_theta_ctx_exp_tau_div_4_inv(ctx), j, k), prec);
                    }
                }
                acb_set(&acb_theta_ctx_exp_tau_a_div_2_inv(ctx, a)[j], x);
                acb_sqr(&acb_theta_ctx_exp_tau_a_inv(ctx, a)[j], x, prec);
            }
            acb_one(x);
            for (j = 0; j < g; j++)
            {
                if (acb_theta_aj_is_zero(a, j, g))
                {
                    continue;
                }
                for (k = j; k < g; k++)
                {
                    if (acb_theta_aj_is_zero(a, k, g))
                    {
                        continue;
                    }
                    acb_mul(x, x, acb_mat_entry(acb_theta_ctx_exp_tau_div_4(ctx), j, k), prec);
                }
            }
            acb_set(acb_theta_ctx_exp_a_tau_a_div_4(ctx, a), x);
        }
    }

    acb_clear(x);
}
