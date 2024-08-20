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
#include "acb_theta.h"

void
acb_theta_ctx_set_z_ql(acb_theta_ctx_t ctx, acb_srcptr z, slong prec)
{
    slong g = acb_theta_ctx_g(ctx);
    slong n = 1 << g;
    slong lp = ACB_THETA_LOW_PREC;
    int z_is_zero = _acb_vec_is_zero(z, g);
    int z_is_real = _acb_vec_is_real(z, g);
    acb_ptr zero;
    arb_ptr w;
    slong a;

    zero = _acb_vec_init(g);
    w = _arb_vec_init(g);

    /* Set exponentials, etc. */
    acb_theta_ctx_set_z(ctx, zero, 0, prec);
    if (!z_is_zero)
    {
        acb_theta_ctx_set_z(ctx, z, 3, prec);
    }

    /* Compute distances */
    if (g >= 2)
    {
        for (a = 0; a < n; a++)
        {
            acb_theta_char_get_arb(w, a, g);
            arb_mat_vector_mul_col(w, acb_theta_ctx_cho(ctx), w, lp);
            acb_theta_dist_lat(&acb_theta_ctx_d0(ctx)[a], w, acb_theta_ctx_cho(ctx), lp);
        }
        if (!z_is_real)
        {
            for (a = 0; a < n; a++)
            {
                acb_theta_char_get_arb(w, a, g);
                arb_mat_vector_mul_col(w, acb_theta_ctx_cho(ctx), w, lp);
                _arb_vec_add(w, acb_theta_ctx_vs(ctx) + 3 * g, w, g, lp);
                acb_theta_dist_lat(&acb_theta_ctx_d(ctx)[a], w, acb_theta_ctx_cho(ctx), lp);
            }
        }
    }

    /* Set info */
    ctx->t_is_zero = 1;
    ctx->z_is_zero = z_is_zero;
    ctx->z_is_real = z_is_real;

    _acb_vec_clear(zero, g);
    _arb_vec_clear(w, g);
}
