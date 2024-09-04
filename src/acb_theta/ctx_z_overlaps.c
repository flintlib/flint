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

    res =  _acb_vec_overlaps(ctx1->exp_z, ctx2->exp_z, g)
        && (ctx1->is_real == ctx2->is_real)
        && arb_overlaps(&ctx1->u, &ctx2->u)
        && arb_overlaps(&ctx1->uinv, &ctx2->uinv);

    if (g > 1)
    {
        res = res && _acb_vec_overlaps(ctx1->exp_2z, ctx2->exp_2z, g)
            && _acb_vec_overlaps(ctx1->exp_z_inv, ctx2->exp_z_inv, g)
            && _acb_vec_overlaps(ctx1->exp_2z_inv, ctx2->exp_2z_inv, g)
            && _arb_vec_overlaps(ctx1->v, ctx2->v, g);
    }

    return res;
}
