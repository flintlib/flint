/*
    Copyright (C) 2025 Jean Kieffer

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

int
acb_theta_ctx_tau_overlaps(const acb_theta_ctx_tau_t ctx1, const acb_theta_ctx_tau_t ctx2)
{
    slong g = ctx2->g;
    slong n = 1 << g;
    int res;

    FLINT_ASSERT(ctx1->g == g);
    FLINT_ASSERT(ctx1->allow_shift == ctx2->allow_shift);

    res = arb_mat_overlaps(&ctx1->yinv, &ctx2->yinv)
        && arb_mat_overlaps(&ctx1->cho, &ctx2->cho)
        && acb_mat_overlaps(ctx1->exp_tau_div_4, ctx2->exp_tau_div_4)
        && acb_mat_overlaps(ctx1->exp_tau_div_2, ctx2->exp_tau_div_2)
        && acb_mat_overlaps(ctx1->exp_tau, ctx2->exp_tau)
        && acb_mat_overlaps(ctx1->exp_tau_div_4_inv, ctx2->exp_tau_div_4_inv)
        && acb_mat_overlaps(ctx1->exp_tau_div_2_inv, ctx2->exp_tau_div_2_inv)
        && acb_mat_overlaps(ctx1->exp_tau_inv, ctx2->exp_tau_inv);

    if (ctx1->allow_shift && res)
    {
        res = _acb_vec_overlaps(ctx1->exp_tau_a, ctx2->exp_tau_a, n * g)
            && _acb_vec_overlaps(ctx1->exp_tau_a_inv, ctx2->exp_tau_a_inv, n * g)
            && _acb_vec_overlaps(ctx1->exp_a_tau_a_div_4, ctx2->exp_a_tau_a_div_4, n);
    }

    return res;
}
