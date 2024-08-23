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
acb_theta_ctx_z_add_real(acb_theta_ctx_z_t res, const acb_theta_ctx_z_t ctx1,
    const acb_theta_ctx_z_t ctx_real, slong prec)
{
    slong g = ctx1->g;
    arb_ptr t_real;
    acb_t x;
    slong j;

    FLINT_ASSERT(ctx_real->g == g);
    FLINT_ASSERT(res->g == g);

    t_real = _arb_vec_init(g);
    acb_init(x);

    /* Copy things */
    _acb_vec_get_real(t_real, acb_theta_ctx_z(ctx_real), g);
    _acb_vec_add(acb_theta_ctx_z(res), acb_theta_ctx_z(ctx1),
        acb_theta_ctx_z(ctx_real), g, prec);
    _arb_vec_set(acb_theta_ctx_r(res), acb_theta_ctx_r(ctx1), g);
    acb_theta_ctx_is_real(res) = acb_theta_ctx_is_real(ctx1);
    if (g > 1)
    {
        _arb_vec_set(acb_theta_ctx_v(res), acb_theta_ctx_v(ctx1), g);
        arb_set(acb_theta_ctx_u(res), acb_theta_ctx_u(ctx1));
    }

    /* Exponentials */
    for (j = 0; j < g; j++)
    {
        acb_mul(&acb_theta_ctx_exp_z(res)[j], &acb_theta_ctx_exp_z(ctx1)[j],
            &acb_theta_ctx_exp_z(ctx_real)[j], prec);
        if (g > 1)
        {
            acb_mul(&acb_theta_ctx_exp_z_inv(res)[j], &acb_theta_ctx_exp_z_inv(ctx1)[j],
                &acb_theta_ctx_exp_z_inv(ctx_real)[j], prec);
            acb_sqr(&acb_theta_ctx_exp_2z(res)[j],
                &acb_theta_ctx_exp_z(res)[j], prec);
            acb_theta_ctx_sqr_inv(&acb_theta_ctx_exp_2z_inv(res)[j],
                &acb_theta_ctx_exp_z_inv(res)[j], &acb_theta_ctx_exp_2z(res)[j],
                acb_theta_ctx_is_real(ctx1), prec);
        }
    }

    /* The factor c gets multiplied by exp(-2 pi i r^T t) */
    arb_dot(acb_realref(x), NULL, 1, acb_theta_ctx_r(ctx1), 1, t_real, 1, g, prec);
    acb_mul_2exp_si(x, x, 1);
    acb_exp_pi_i(x, x, prec);
    acb_mul(acb_theta_ctx_c(res), acb_theta_ctx_c(ctx1), x, prec);

    acb_clear(x);
    _arb_vec_clear(t_real, g);
}
