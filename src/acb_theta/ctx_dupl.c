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
acb_theta_ctx_dupl(acb_theta_ctx_t ctx, slong prec)
{
    slong g = acb_theta_ctx_g(ctx);
    slong n = 1 << g;
    slong nb = (ctx->t_is_zero ? 1 : 3);
    arb_t sqrt2;
    slong j, k;

    arb_init(sqrt2);
    arb_set_si(sqrt2, 2);
    arb_sqrt(sqrt2, sqrt2, prec);

    /* Duplicate tau part */
    acb_mat_scalar_mul_2exp_si(acb_theta_ctx_tau(ctx), acb_theta_ctx_tau(ctx), 1);
    arb_mat_scalar_mul_2exp_si(acb_theta_ctx_y(ctx), acb_theta_ctx_y(ctx), 1);
    arb_mat_scalar_mul_2exp_si(acb_theta_ctx_yinv(ctx), acb_theta_ctx_yinv(ctx), -1);
    acb_mat_set(acb_theta_ctx_exp_tau_div_4(ctx), acb_theta_ctx_exp_tau_div_2(ctx));
    acb_mat_set(acb_theta_ctx_exp_tau_div_2(ctx), acb_theta_ctx_exp_tau(ctx));
    for (j = 0; j < g; j++)
    {
        for (k = j; k < g; k++)
        {
            acb_sqr(acb_mat_entry(acb_theta_ctx_exp_tau(ctx), j, k),
                acb_mat_entry(acb_theta_ctx_exp_tau(ctx), j, k), prec);
        }
    }
    if (g > 1)
    {
        arb_mat_scalar_mul_arb(acb_theta_ctx_cho(ctx), acb_theta_ctx_cho(ctx), sqrt2, prec);
        arb_mat_scalar_div_arb(acb_theta_ctx_choinv(ctx), acb_theta_ctx_choinv(ctx), sqrt2, prec);
        for (j = 0; j < g; j++)
        {
            for (k = j + 1; k < g; k++)
            {
                acb_sqr(acb_mat_entry(acb_theta_ctx_exp_tau_inv(ctx), j, k),
                    acb_mat_entry(acb_theta_ctx_exp_tau_inv(ctx), j, k), prec);
            }
        }
    }

    /* Duplicate z part. Do nothing for the zero vector. */
    /* For the completely real part, c, v, a, u do not change */
    for (j = 1; j < nb; j++)
    {
        _acb_vec_set(acb_theta_ctx_exp_zs(ctx) + j * g, acb_theta_ctx_exp_2zs(ctx) + j * g, g);
        _acb_vec_set(acb_theta_ctx_exp_zs_inv(ctx) + j * g, acb_theta_ctx_exp_2zs_inv(ctx) + j * g, g);
    }
    if (!ctx->t_is_zero)
    {
        _acb_vec_set(acb_theta_ctx_exp_2zs(ctx) + g, acb_theta_ctx_exp_2zs(ctx) + 2 * g, g);
        _acb_vec_set(acb_theta_ctx_exp_2zs_inv(ctx) + g, acb_theta_ctx_exp_2zs_inv(ctx) + 2 * g, g);
        for (j = 0; j < g; j++)
        {
            acb_sqr(&acb_theta_ctx_exp_2zs(ctx)[2 * g + j], &acb_theta_ctx_exp_2zs(ctx)[2 * g + j], prec);
            acb_conj(&acb_theta_ctx_exp_2zs_inv(ctx)[2 * g + j], &acb_theta_ctx_exp_2zs(ctx)[2 * g + j]);
        }
    }
    /* For the non-real part, c, v, and u may change */
    if (!ctx->z_is_zero)
    {
        /* Set exponentials and cs */
        for (j = 3; j < 3 + nb; j++)
        {
            _acb_vec_set(acb_theta_ctx_exp_zs(ctx) + j * g, acb_theta_ctx_exp_2zs(ctx) + j * g, g);
            _acb_vec_set(acb_theta_ctx_exp_zs_inv(ctx) + j * g, acb_theta_ctx_exp_2zs_inv(ctx) + j * g, g);
            acb_sqr(&acb_theta_ctx_cs(ctx)[j], &acb_theta_ctx_cs(ctx)[j], prec);
        }
        for (j = 0; j < nb * g; j++)
        {
            acb_sqr(&acb_theta_ctx_exp_2zs(ctx)[3 * g + j], &acb_theta_ctx_exp_2zs(ctx)[3 * g + j], prec);
            if (ctx->z_is_real)
            {
                acb_conj(&acb_theta_ctx_exp_2zs_inv(ctx)[3 * g + j], &acb_theta_ctx_exp_2zs(ctx)[3 * g + j]);
            }
            else
            {
                acb_sqr(&acb_theta_ctx_exp_2zs_inv(ctx)[3 * g + j],
                    &acb_theta_ctx_exp_2zs_inv(ctx)[3 * g + j], prec);
            }
        }

        /* Set us and vs (same for all 3 vectors) */
        arb_sqr(&acb_theta_ctx_us(ctx)[3], &acb_theta_ctx_us(ctx)[3], prec);
        if (g > 1)
        {
            _arb_vec_scalar_mul(acb_theta_ctx_vs(ctx) + 3 * g, acb_theta_ctx_vs(ctx) + 3 * g, g, sqrt2, prec);
        }
        for (j = 4; j < 3 + nb; j++)
        {
            arb_set(&acb_theta_ctx_us(ctx)[j], &acb_theta_ctx_us(ctx)[3]);
            if (g > 1)
            {
                _arb_vec_set(acb_theta_ctx_vs(ctx) + j * g, acb_theta_ctx_vs(ctx) + 3 * g, g);
            }
        }
    }

    /* Duplicate distances */
    if (g > 1)
    {
        _arb_vec_scalar_mul_2exp_si(acb_theta_ctx_d0(ctx), acb_theta_ctx_d0(ctx), n, 1);
        if (!ctx->z_is_real)
        {
            _arb_vec_scalar_mul_2exp_si(acb_theta_ctx_d(ctx), acb_theta_ctx_d(ctx), n, 1);
        }
    }

    arb_clear(sqrt2);
}
