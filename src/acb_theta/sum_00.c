/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "acb.h"
#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_modular.h"
#include "acb_theta.h"

static void
acb_theta_sum_00_worker(acb_ptr th, acb_srcptr v1, acb_srcptr v2,
    const slong * precs, slong len, const acb_t cofactor, const slong * coords,
    slong ord, slong g, slong prec, slong fullprec)
{
    acb_t sum;

    acb_init(sum);

    acb_dot(sum, NULL, 0, v1, 1, v2, 1, len, prec);
    acb_addmul(th, sum, cofactor, fullprec);

    acb_clear(sum);
}

void
acb_theta_sum_00(acb_ptr th, const acb_theta_ctx_z_struct * vec, slong nb,
    const acb_theta_ctx_tau_t ctx_tau, slong prec)
{
    slong g = ctx_tau->g;
    slong j;

    FLINT_ASSERT(nb >= 0);
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
            acb_modular_theta_sum(&res[0], &res[1], &res[2], &res[3],
                (&vec[j])->exp_z, (&vec[j])->is_real,
                acb_mat_entry(ctx_tau->exp_tau, 0, 0), 1, prec);
            acb_set(&th[j], &res[2]);
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

        acb_theta_ctx_z_common_v(v, vec, nb, prec);
        acb_theta_sum_radius(R2, eps, &ctx_tau->cho, 0, prec);
        b = acb_theta_eld_set(E, &ctx_tau->cho, R2, v);

        if (b)
        {
            for (j = 0; j < nb; j++)
            {
                acb_theta_sum_work(&th[j], 1, (&vec[j])->exp_2z,
                    (&vec[j])->exp_2z_inv, 1,
                    ctx_tau->exp_tau, ctx_tau->exp_tau_inv, E, 0,
                    prec, acb_theta_sum_00_worker);
                arb_mul_arf(err, &(&vec[j])->u, eps, prec);
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
