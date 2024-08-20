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
#include "acb_modular.h"
#include "acb_theta.h"

void
acb_theta_sum_00(acb_ptr th, const acb_theta_ctx_t ctx, slong prec)
{
    slong g = acb_theta_ctx_g(ctx);
    slong nb = acb_theta_ctx_nb(ctx);
    slong j;

    if (nb == 0)
    {
        return;
    }

    if (g == 1)
    {
        acb_ptr res;

        res = _acb_vec_init(4);
        for (j = 0; j < nb; j++)
        {
            /* acb_modular_theta_sum recomputes the inverse of exp_z */
            /* todo: store w_is_unit as part of context */
            acb_modular_theta_sum(&res[0], &res[1], &res[2], &res[3],
                &acb_theta_ctx_exp_zs(ctx)[j], 0,
                acb_mat_entry(acb_theta_ctx_exp_tau(ctx), 0, 0), 1, prec);
            acb_mul(&th[j], &res[2], &acb_theta_ctx_cs(ctx)[j], prec);
        }
        _acb_vec_clear(res, 4);
    }
    else
    {
        acb_theta_eld_t E;
        arf_t R2, eps;
        arb_t err;
        arb_ptr v;
        int b;

        acb_theta_eld_init(E, g, g);
        arf_init(R2);
        arf_init(eps);
        arb_init(err);
        v = _arb_vec_init(g);

        acb_theta_ctx_common_v(v, ctx, prec);
        acb_theta_naive_radius(R2, eps, acb_theta_ctx_cho(ctx), 0, prec);
        b = acb_theta_eld_set(E, acb_theta_ctx_cho(ctx), R2, v);

        if (b)
        {
            acb_theta_sum_work(th, 1, acb_theta_ctx_exp_2zs(ctx), acb_theta_ctx_exp_2zs_inv(ctx), nb,
                acb_theta_ctx_exp_tau(ctx), acb_theta_ctx_exp_tau_inv(ctx), E, 0,
                prec, acb_theta_sum_00_worker);
            for (j = 0; j < nb; j++)
            {
                acb_mul(&th[j], &th[j], &acb_theta_ctx_cs(ctx)[j], prec);
                arb_mul_arf(err, &acb_theta_ctx_us(ctx)[j], eps, prec);
                acb_add_error_arb(&th[j], err);
            }
        }
        else
        {
            for (j = 0; j < nb; j++)
            {
                acb_indeterminate(&th[j]);
            }
        }

        acb_theta_eld_clear(E);
        arf_clear(R2);
        arf_clear(eps);
        arb_clear(err);
        _arb_vec_clear(v, g);
    }
}
