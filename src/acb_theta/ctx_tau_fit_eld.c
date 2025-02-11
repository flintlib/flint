/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"

void
acb_theta_ctx_tau_fit_eld(acb_theta_ctx_tau_t ctx, const acb_theta_eld_t E, slong prec)
{
    slong g = ctx->g;
    acb_ptr cpy;
    acb_t d, dd;
    slong j, k, len;

    acb_init(d);
    acb_init(dd);

    for (j = 0; j < g; j++)
    {
        len = ctx->sqr_pow_len[j];
        if (len < (E->box[j] + 1))
        {
            cpy = ctx->sqr_pow[j];
            ctx->sqr_pow[j] = _acb_vec_init(E->box[j] + 1);
            ctx->sqr_pow_len[j] = (E->box[j] + 1);
            _acb_vec_set(ctx->sqr_pow[j], cpy, len);
            _acb_vec_clear(cpy, len);

            acb_sqr(dd, acb_mat_entry(ctx->exp_tau, j, j), prec);
            if (len == 0)
            {
                acb_one(&(ctx->sqr_pow[j])[len]);
                acb_set(d, acb_mat_entry(ctx->exp_tau, j, j));
            }
            else
            {
                acb_pow_ui(d, acb_mat_entry(ctx->exp_tau, j, j), 2 * len - 1, prec);
                acb_mul(&(ctx->sqr_pow[j])[len], &(ctx->sqr_pow[j])[len - 1], d, prec);
                acb_mul(d, d, dd, prec);
            }

            for (k = len + 1; k <= (E->box[j]); k++)
            {
                acb_mul(&(ctx->sqr_pow[j])[k], &(ctx->sqr_pow[j])[k - 1], d, prec);
                acb_mul(d, d, dd, prec);
            }
        }
    }

    acb_clear(d);
    acb_clear(dd);
}
