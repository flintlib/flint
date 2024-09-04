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

static void
acb_theta_sum_0b_worker(acb_ptr th, acb_srcptr v1, acb_srcptr v2,
    const slong * precs, slong len, const acb_t cofactor, const slong * coords,
    slong ord, slong g, slong prec, slong fullprec)
{
    slong n = 1 << g;
    acb_t s0, s1, add, sub;
    ulong b;
    slong dot;

    acb_init(s0);
    acb_init(s1);
    acb_init(add);
    acb_init(sub);

    /* Compute alternate sums to adjust signs */
    acb_dot(s0, NULL, 0, v1, 2, v2, 2, (len + 1) / 2, prec);
    acb_dot(s1, NULL, 0, v1 + 1, 2, v2 + 1, 2, len / 2, prec);
    acb_add(add, s0, s1, prec);
    acb_sub(sub, s0, s1, prec);
    acb_mul(add, add, cofactor, prec);
    acb_mul(sub, sub, cofactor, prec);

    for (b = 0; b < n; b++)
    {
        dot = acb_theta_char_dot_slong(b, coords, g) % 2;
        if ((b >> (g - 1)) && dot)
        {
            acb_sub(&th[b], &th[b], sub, fullprec);
        }
        else if ((b >> (g - 1)))
        {
            acb_add(&th[b], &th[b], sub, fullprec);
        }
        else if (dot)
        {
            acb_sub(&th[b], &th[b], add, fullprec);
        }
        else
        {
            acb_add(&th[b], &th[b], add, fullprec);
        }
    }

    acb_clear(s0);
    acb_clear(s1);
    acb_clear(add);
    acb_clear(sub);
}

void
acb_theta_sum_0b(acb_ptr th, const acb_theta_ctx_z_struct * vec, slong nb,
    const acb_theta_ctx_tau_t ctx_tau, slong prec)
{
    slong g = ctx_tau->g;
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
                (&vec[j])->exp_z, (&vec[j])->is_real,
                acb_mat_entry(ctx_tau->exp_tau, 0, 0), 1, prec);
            _acb_vec_set(th + 2 * j, res + 2, 2);
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
                acb_theta_sum_work(th + j * n, n, (&vec[j])->exp_2z,
                    (&vec[j])->exp_2z_inv, 1,
                    ctx_tau->exp_tau, ctx_tau->exp_tau_inv, E, 0,
                    prec, acb_theta_sum_0b_worker);
                arb_mul_arf(err, &(&vec[j])->u, eps, prec);
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
