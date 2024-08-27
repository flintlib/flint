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

int
acb_theta_ctx_tau_overlaps(const acb_theta_ctx_tau_t ctx1, const acb_theta_ctx_tau_t ctx2)
{
    slong g = acb_theta_ctx_g(ctx2);
    slong n = 1 << g;
    int res;

    FLINT_ASSERT(acb_theta_ctx_g(ctx1) == g);

    res = acb_mat_overlaps(acb_theta_ctx_tau(ctx1), acb_theta_ctx_tau(ctx2))
        && arb_mat_overlaps(acb_theta_ctx_y(ctx1),
            acb_theta_ctx_y(ctx2))
        && arb_mat_overlaps(acb_theta_ctx_yinv(ctx1),
            acb_theta_ctx_yinv(ctx2))
        && acb_mat_overlaps(acb_theta_ctx_exp_tau_div_4(ctx1),
            acb_theta_ctx_exp_tau_div_4(ctx2))
        && acb_mat_overlaps(acb_theta_ctx_exp_tau_div_2(ctx1),
            acb_theta_ctx_exp_tau_div_2(ctx2))
        && acb_mat_overlaps(acb_theta_ctx_exp_tau(ctx1),
            acb_theta_ctx_exp_tau(ctx2));

    if (g > 1 && res)
    {
        res = arb_mat_overlaps(acb_theta_ctx_cho(ctx1), acb_theta_ctx_cho(ctx2))
            && arb_mat_overlaps(acb_theta_ctx_choinv(ctx1),
                acb_theta_ctx_choinv(ctx2))
            && acb_mat_overlaps(acb_theta_ctx_exp_tau_div_4_inv(ctx1),
                acb_theta_ctx_exp_tau_div_4_inv(ctx2))
            && acb_mat_overlaps(acb_theta_ctx_exp_tau_div_2_inv(ctx1),
                acb_theta_ctx_exp_tau_div_2_inv(ctx2))
            && acb_mat_overlaps(acb_theta_ctx_exp_tau_inv(ctx1),
                acb_theta_ctx_exp_tau_inv(ctx2))
            && _acb_vec_overlaps(ctx1->exp_tau_a_div_2,
                ctx2->exp_tau_a_div_2, n * g)
            && _acb_vec_overlaps(ctx1->exp_tau_a,
                ctx2->exp_tau_a, n * g)
            && _acb_vec_overlaps(ctx1->exp_tau_a_div_2_inv,
                ctx2->exp_tau_a_div_2_inv, n * g)
            && _acb_vec_overlaps(ctx1->exp_tau_a_inv,
                ctx2->exp_tau_a_inv, n * g)
            && _acb_vec_overlaps(ctx1->exp_a_tau_a_div_4,
                ctx2->exp_a_tau_a_div_4, n);
    }

    return res;
}
