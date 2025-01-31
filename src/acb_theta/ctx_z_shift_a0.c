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
acb_theta_ctx_z_shift_a0(acb_theta_ctx_z_t res, acb_t c, const acb_theta_ctx_z_t ctx,
    const acb_theta_ctx_tau_t ctx_tau, ulong a, slong prec)
{
    slong g = ctx_tau->g;
    arb_ptr v_shift;
    acb_t cinv;
    arb_t abs;
    slong j;

    v_shift = _arb_vec_init(g);
    acb_init(cinv);
    arb_init(abs);

    /* Do not set exp_z or exp_z_inv. */
    /* Replace exp_2z by analogs for z + tau a/2 */
    for (j = 0; j < g; j++)
    {
        acb_mul(&res->exp_2z[j], &ctx->exp_2z[j],
            &ctx_tau->exp_tau_a[a * g + j], prec);
        acb_mul(&res->exp_2z_inv[j], &ctx->exp_2z_inv[j],
            &ctx_tau->exp_tau_a_inv[a * g + j], prec);
    }

    /* Compute cofactor exp(pi i a^T z), and multiply by common cofactor
       exp(pi i/4 a^T tau a) */
    acb_one(c);
    for (j = 0; j < g; j++)
    {
        if (acb_theta_aj_is_zero(a, j, g))
        {
            continue;
        }
        acb_mul(c, c, &ctx->exp_z[j], prec);
    }
    acb_mul(c, c, &ctx_tau->exp_a_tau_a_div_4[a], prec);

    /* Compute v; u and uinv must be multiplied by abs(c) */
    acb_abs(abs, c, prec);
    arb_mul(&res->uinv, &ctx->uinv, abs, prec);

    arb_inv(abs, abs, prec);
    if (acb_is_finite(c) && !arb_is_finite(abs))
    {
        /* Recompute cinv by multiplications */
        acb_one(cinv);
        for (j = 0; j < g; j++)
        {
            if (acb_theta_aj_is_zero(a, j, g))
            {
                continue;
            }
            acb_mul(cinv, cinv, &ctx->exp_z_inv[j], prec);
        }
        acb_div(cinv, cinv, &ctx_tau->exp_a_tau_a_div_4[a], prec);
        acb_abs(abs, cinv, prec);
    }
    arb_mul(&res->u, &ctx->u, abs, prec);

    acb_theta_char_get_arb(v_shift, a, g);
    arb_mat_vector_mul_col(v_shift, &ctx_tau->cho, v_shift, prec);
    _arb_vec_add(res->v, v_shift, ctx->v, g, prec);

    _arb_vec_clear(v_shift, g);
    acb_clear(cinv);
    arb_clear(abs);
}
