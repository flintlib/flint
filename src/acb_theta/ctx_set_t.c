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
acb_theta_ctx_set_t(acb_theta_ctx_t ctx, const acb_ptr t, slong prec)
{
    slong g = acb_theta_ctx_g(ctx);
    acb_t x;
    arb_ptr t_real;
    slong j;

    FLINT_ASSERT(_acb_vec_is_real(t, g));
    ctx->t_is_zero = _acb_vec_is_zero(t, g);
    if (ctx->t_is_zero)
    {
        return;
    }

    acb_init(x);
    t_real = _arb_vec_init(g);
    acb_theta_ctx_set_z(ctx, t, 1, prec);

    /* Propagate 2t using squares/conjs */
    _acb_vec_set(acb_theta_ctx_exp_zs(ctx) + 2 * g, acb_theta_ctx_exp_2zs(ctx) + g, g);
    _acb_vec_set(acb_theta_ctx_exp_zs_inv(ctx) + 2 * g, acb_theta_ctx_exp_2zs_inv(ctx) + g, g);
    for (j = 0; j < g; j++)
    {
        acb_sqr(&acb_theta_ctx_exp_2zs(ctx)[2 * g + j],
            &acb_theta_ctx_exp_2zs(ctx)[g + j], prec);
        acb_conj(&acb_theta_ctx_exp_2zs_inv(ctx)[2 * g + j],
            &acb_theta_ctx_exp_2zs(ctx)[2 * g + j]);
        acb_one(&acb_theta_ctx_cs(ctx)[2]);
        arb_one(&acb_theta_ctx_us(ctx)[2]);
    }

    /* Propagate z+t, z+2t using mults */
    if (!ctx->z_is_zero)
    {
        for (j = 0; j < g; j++)
        {
            acb_mul(&acb_theta_ctx_exp_zs(ctx)[4 * g + j],
                &acb_theta_ctx_exp_zs(ctx)[3 * g + j],
                &acb_theta_ctx_exp_zs(ctx)[g + j], prec);
            acb_mul(&acb_theta_ctx_exp_zs(ctx)[5 * g + j],
                &acb_theta_ctx_exp_zs(ctx)[3 * g + j],
                &acb_theta_ctx_exp_zs(ctx)[2 * g + j], prec);
            acb_sqr(&acb_theta_ctx_exp_2zs(ctx)[4 * g + j],
                &acb_theta_ctx_exp_zs(ctx)[4 * g + j], prec);
            acb_sqr(&acb_theta_ctx_exp_2zs(ctx)[5 * g + j],
                &acb_theta_ctx_exp_zs(ctx)[5 * g + j], prec);
            if (ctx->z_is_real)
            {
                acb_conj(&acb_theta_ctx_exp_zs_inv(ctx)[4 * g + j],
                    &acb_theta_ctx_exp_zs(ctx)[4 * g + j]);
                acb_conj(&acb_theta_ctx_exp_zs_inv(ctx)[5 * g + j],
                    &acb_theta_ctx_exp_zs(ctx)[5 * g + j]);
                acb_conj(&acb_theta_ctx_exp_2zs_inv(ctx)[4 * g + j],
                    &acb_theta_ctx_exp_2zs(ctx)[4 * g + j]);
                acb_conj(&acb_theta_ctx_exp_2zs_inv(ctx)[5 * g + j],
                    &acb_theta_ctx_exp_2zs(ctx)[5 * g + j]);
            }
            else
            {
                acb_mul(&acb_theta_ctx_exp_zs_inv(ctx)[4 * g + j],
                    &acb_theta_ctx_exp_zs_inv(ctx)[3 * g + j],
                    &acb_theta_ctx_exp_zs_inv(ctx)[g + j], prec);
                acb_mul(&acb_theta_ctx_exp_zs_inv(ctx)[5 * g + j],
                    &acb_theta_ctx_exp_zs_inv(ctx)[3 * g + j],
                    &acb_theta_ctx_exp_zs_inv(ctx)[2 * g + j], prec);
                acb_sqr(&acb_theta_ctx_exp_2zs_inv(ctx)[4 * g + j],
                    &acb_theta_ctx_exp_zs_inv(ctx)[4 * g + j], prec);
                acb_sqr(&acb_theta_ctx_exp_2zs_inv(ctx)[5 * g + j],
                    &acb_theta_ctx_exp_zs_inv(ctx)[5 * g + j], prec);
            }
        }
        /* The factor c gets multiplied by exp(-2 pi i a^T t) */
        _acb_vec_get_real(t_real, t, g);
        arb_dot(acb_realref(x), NULL, 1, acb_theta_ctx_as(ctx) + 3 * g, 1, t_real, 1, g, prec);
        acb_mul_2exp_si(x, x, 1);
        acb_exp_pi_i(x, x, prec);
        acb_mul(&acb_theta_ctx_cs(ctx)[4], &acb_theta_ctx_cs(ctx)[3], x, prec);
        acb_mul(&acb_theta_ctx_cs(ctx)[5], &acb_theta_ctx_cs(ctx)[4], x, prec);
        for (j = 4; j < 6; j++)
        {
            arb_set(&acb_theta_ctx_us(ctx)[j], &acb_theta_ctx_us(ctx)[3]);
            _arb_vec_set(acb_theta_ctx_as(ctx) + j * g, acb_theta_ctx_as(ctx) + 3 * g, g);
            if (g > 1)
            {
                _arb_vec_set(acb_theta_ctx_vs(ctx) + j * g, acb_theta_ctx_vs(ctx) + 3 * g, g);
            }
        }
    }

    acb_clear(x);
    _arb_vec_clear(t_real, g);
}
