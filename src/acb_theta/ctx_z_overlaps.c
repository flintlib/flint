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

int acb_theta_ctx_z_overlaps(const acb_theta_ctx_z_t ctx1, const acb_theta_ctx_z_t ctx2)
{
    slong g = ctx1->g;
    int res;

    if (ctx2->g != g)
    {
        return 0;
    }

    res = _acb_vec_overlaps(acb_theta_ctx_z(ctx1),
        acb_theta_ctx_z(ctx2), g)
        && _acb_vec_overlaps(acb_theta_ctx_exp_z(ctx1),
            acb_theta_ctx_exp_z(ctx2), g)
        && acb_overlaps(acb_theta_ctx_c(ctx1),
            acb_theta_ctx_c(ctx2))
        && _arb_vec_overlaps(acb_theta_ctx_r(ctx1),
            acb_theta_ctx_r(ctx2), g)
        && (acb_theta_ctx_is_real(ctx1) == acb_theta_ctx_is_real(ctx2))
        && arb_overlaps(acb_theta_ctx_u(ctx1), acb_theta_ctx_u(ctx2))
        && arb_overlaps(acb_theta_ctx_uinv(ctx1), acb_theta_ctx_uinv(ctx2));

    if (g > 1)
    {
        res = res && _acb_vec_overlaps(acb_theta_ctx_exp_2z(ctx1),
            acb_theta_ctx_exp_2z(ctx2), g)
            && _acb_vec_overlaps(acb_theta_ctx_exp_z_inv(ctx1),
                acb_theta_ctx_exp_z_inv(ctx2), g)
            && _acb_vec_overlaps(acb_theta_ctx_exp_2z_inv(ctx1),
                acb_theta_ctx_exp_2z_inv(ctx2), g)
            && _arb_vec_overlaps(acb_theta_ctx_v(ctx1),
                acb_theta_ctx_v(ctx2), g);
    }

    return res;
}
