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
acb_theta_ctx_shift_z(acb_theta_ctx_t new_ctx, const acb_theta_ctx_t ctx,
    slong start, slong nb, ulong a, slong prec)
{
    slong g = acb_theta_ctx_g(ctx);
    acb_t c, cinv, csqr, csqrinv;
    arb_t abs;
    slong j, k;

    acb_init(c);
    acb_init(cinv);
    acb_init(csqr);
    acb_init(csqrinv);
    arb_init(abs);

    acb_theta_ctx_copy_tau(new_ctx, ctx);

    /* Replace exp_zs by analogs for z + tau a/2 */
    _acb_vec_set(acb_theta_ctx_exp_zs(new_ctx), acb_theta_ctx_exp_zs(ctx) + start * g, nb * g);
    _acb_vec_set(acb_theta_ctx_exp_zs_inv(new_ctx), acb_theta_ctx_exp_zs_inv(ctx) + start * g, nb * g);
    _acb_vec_set(acb_theta_ctx_exp_2zs(new_ctx), acb_theta_ctx_exp_2zs(ctx) + start * g, nb * g);
    _acb_vec_set(acb_theta_ctx_exp_2zs_inv(new_ctx), acb_theta_ctx_exp_2zs_inv(ctx) + start * g, nb * g);
    for (j = 0; j < g; j++)
    {
        acb_one(c);
        for (k = 0; k < g; k++)
        {
            if (acb_theta_aj_is_zero(a, k, g))
            {
                continue;
            }
            if (k < j)
            {
                acb_mul(c, c, acb_mat_entry(acb_theta_ctx_exp_tau_div_4(ctx), k, j), prec);
            }
            else if (k == j)
            {
                acb_mul(c, c, acb_mat_entry(acb_theta_ctx_exp_tau_div_2(ctx), k, k), prec);
            }
            else
            {
                acb_mul(c, c, acb_mat_entry(acb_theta_ctx_exp_tau_div_4(ctx), j, k), prec);
            }
        }
        acb_inv(cinv, c, prec);
        acb_sqr(csqr, c, prec);
        acb_sqr(csqrinv, cinv, prec);
        for (k = 0; k < nb; k++)
        {
            acb_mul(&acb_theta_ctx_exp_zs(new_ctx)[k * g + j],
                &acb_theta_ctx_exp_zs(new_ctx)[k * g + j], c, prec);
            acb_mul(&acb_theta_ctx_exp_zs_inv(new_ctx)[k * g + j],
                &acb_theta_ctx_exp_zs_inv(new_ctx)[k * g + j], cinv, prec);
            acb_mul(&acb_theta_ctx_exp_2zs(new_ctx)[k * g + j],
                &acb_theta_ctx_exp_2zs(new_ctx)[k * g + j], csqr, prec);
            acb_mul(&acb_theta_ctx_exp_2zs_inv(new_ctx)[k * g + j],
                &acb_theta_ctx_exp_2zs_inv(new_ctx)[k * g + j], csqrinv, prec);
        }
    }

    /* For each z, compute cofactor exp(pi i a^T z) */
    for (j = 0; j < nb; j++)
    {
        acb_one(c);
        for (k = 0; k < g; k++)
        {
            if (acb_theta_aj_is_zero(a, k, g))
            {
                continue;
            }
            acb_mul(c, c, &acb_theta_ctx_exp_zs(ctx)[(start + j) * g + k], prec);
        }
        acb_set(&acb_theta_ctx_cs(new_ctx)[j], c);
    }

    /* A miracle happens: ctx_as does not appear in this function. */

    /* Compute common cofactor exp(pi i/4 a^T tau a) */
    acb_one(c);
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
            acb_mul(c, c, acb_mat_entry(acb_theta_ctx_exp_tau_div_4(ctx), j, k), prec);
        }
    }
    _acb_vec_scalar_mul(acb_theta_ctx_cs(new_ctx), acb_theta_ctx_cs(new_ctx), nb, c, prec);

    /* Compute cs, us, as, vs */
    _arb_vec_set(acb_theta_ctx_as(new_ctx), acb_theta_ctx_as(ctx) + start * g, nb * g);
    for (j = 0; j < nb; j++)
    {
        acb_abs(abs, &acb_theta_ctx_cs(new_ctx)[j], prec);
        arb_mul(&acb_theta_ctx_us(new_ctx)[j], &acb_theta_ctx_us(ctx)[start + j], abs, prec);
        acb_mul(&acb_theta_ctx_cs(new_ctx)[j], &acb_theta_ctx_cs(new_ctx)[j], &acb_theta_ctx_cs(ctx)[start + j], prec);
    }
    if (g > 1)
    {
        arb_ptr v_shift;
        v_shift = _arb_vec_init(g);

        acb_theta_char_get_arb(v_shift, a, g);
        arb_mat_vector_mul_col(v_shift, acb_theta_ctx_cho(ctx), v_shift, prec);
        for (j = 0; j < nb; j++)
        {
            _arb_vec_add(acb_theta_ctx_vs(new_ctx) + j * g, v_shift,
                acb_theta_ctx_vs(ctx) + (start + j) * g, g, prec);
        }

        _arb_vec_clear(v_shift, g);
    }

    acb_clear(c);
    acb_clear(cinv);
    acb_clear(csqr);
    acb_clear(csqrinv);
    arb_clear(abs);
}
