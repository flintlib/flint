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
acb_theta_sum_0b(acb_ptr th, const acb_theta_ctx_z_struct * vec, slong nb,
    const acb_theta_ctx_tau_t ctx_tau, slong prec)
{
    slong g = acb_theta_ctx_g(ctx_tau);
    slong n = 1 << g;
    slong j, k;

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
                acb_theta_ctx_exp_z(&vec[j]), acb_theta_ctx_is_real(&vec[j]),
                acb_mat_entry(acb_theta_ctx_exp_tau(ctx_tau), 0, 0), 1, prec);
            acb_mul(&th[2 * j], &res[2], acb_theta_ctx_c(&vec[j]), prec);
            acb_mul(&th[2 * j + 1], &res[3], acb_theta_ctx_c(&vec[j]), prec);
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
        acb_theta_naive_radius(R2, eps, acb_theta_ctx_cho(ctx_tau), 0, prec);
        b = acb_theta_eld_set(E, acb_theta_ctx_cho(ctx_tau), R2, v);

        if (b)
        {
            for (j = 0; j < nb; j++)
            {
                acb_theta_sum_work(th + j * n, n, acb_theta_ctx_exp_2z(&vec[j]),
                    acb_theta_ctx_exp_2z_inv(&vec[j]), 1,
                    acb_theta_ctx_exp_tau(ctx_tau), acb_theta_ctx_exp_tau_inv(ctx_tau), E, 0,
                    prec, acb_theta_sum_0b_worker);
                _acb_vec_scalar_mul(th + j * n, th + j * n, n, acb_theta_ctx_c(&vec[j]), prec);

                arb_mul_arf(err, acb_theta_ctx_u(&vec[j]), eps, prec);
                for (k = 0; k < n; k++)
                {
                    acb_add_error_arb(&th[j * n + k], err);
                }
            }
        }
        else
        {
            for (j = 0; j < nb * n; j++)
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
