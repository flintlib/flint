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

void
acb_theta_ctx_clear(acb_theta_ctx_t ctx)
{
    slong g = acb_theta_ctx_g(ctx);
    slong nb = acb_theta_ctx_nb(ctx);
    slong n = 1 << g;

    acb_mat_clear(acb_theta_ctx_tau(ctx));
    arb_mat_clear(acb_theta_ctx_y(ctx));
    arb_mat_clear(acb_theta_ctx_yinv(ctx));
    acb_mat_clear(acb_theta_ctx_exp_tau_div_4(ctx));
    acb_mat_clear(acb_theta_ctx_exp_tau_div_2(ctx));
    acb_mat_clear(acb_theta_ctx_exp_tau(ctx));
    _acb_vec_clear(acb_theta_ctx_exp_zs(ctx), nb * g);
    _acb_vec_clear(acb_theta_ctx_exp_zs_inv(ctx), nb * g);
    _acb_vec_clear(acb_theta_ctx_exp_2zs(ctx), nb * g);
    _acb_vec_clear(acb_theta_ctx_exp_2zs_inv(ctx), nb * g);
    _acb_vec_clear(acb_theta_ctx_cs(ctx), nb);
    _arb_vec_clear(acb_theta_ctx_us(ctx), nb);
    _arb_vec_clear(acb_theta_ctx_as(ctx), nb * g);

    if (g >= 2)
    {
        arb_mat_clear(acb_theta_ctx_cho(ctx));
        arb_mat_clear(acb_theta_ctx_choinv(ctx));
        acb_mat_clear(acb_theta_ctx_exp_tau_inv(ctx));
        _arb_vec_clear(acb_theta_ctx_vs(ctx), nb * g);
        _arb_vec_clear(acb_theta_ctx_d0(ctx), n);
        _arb_vec_clear(acb_theta_ctx_d(ctx), n);
    }
}
