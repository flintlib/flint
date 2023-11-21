/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"
#include "acb_modular.h"
#include "acb_theta.h"

static int
acb_theta_ql_a0_naive_gen(acb_ptr th, acb_srcptr t, acb_srcptr z, arb_srcptr d0,
    arb_srcptr d, const acb_mat_t tau, slong guard, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    int hast = !_acb_vec_is_zero(t, g);
    int hasz = !_acb_vec_is_zero(z, g);
    slong nbt = (hast ? 3 : 1);
    slong nbz = (hasz ? 2 : 1);
    acb_ptr x, aux;
    slong j, k;
    int res;

    x = _acb_vec_init(g * nbt);
    aux = _acb_vec_init(nbt);

    for (k = 0; k < nbt; k++)
    {
        _acb_vec_scalar_mul_ui(x + k * g, t, g, k, prec);
    }
    for (k = 0; k < n; k++)
    {
        acb_theta_naive_fixed_ab(aux, k << g, x, nbt, tau,
            prec + acb_theta_dist_addprec(&d0[k]));
        for (j = 0; j < nbt; j++)
        {
            acb_set(&th[j * n + k], &aux[j]);
        }
    }

    if (hasz)
    {
        for (k = 0; k < nbt; k++)
        {
            _acb_vec_add(x + k * g, x + k * g, z, g, prec);
        }
        for (k = 0; k < n; k++)
        {
            acb_theta_naive_fixed_ab(aux, k << g, x, nbt, tau,
                prec + acb_theta_dist_addprec(&d[k]));
            for (j = 0; j < nbt; j++)
            {
                acb_set(&th[nbt * n + j * n + k], &aux[j]);
            }
        }
    }
    res = _acb_vec_is_finite(th, n * nbz * nbt);

    _acb_vec_clear(x, g * nbt);
    _acb_vec_clear(aux, nbt);
    return res;
}

/* when g = 1, we don't go as far near the cusp and computing exponentials is
   relatively more expensive */
static int
acb_theta_ql_a0_naive_g1(acb_ptr th, acb_srcptr t, acb_srcptr z, arb_srcptr d0,
    arb_srcptr d, const acb_mat_t tau, slong guard, slong prec)
{
    int hast = !acb_is_zero(t);
    int hasz = !acb_is_zero(z);
    slong nbt = (hast ? 3 : 1);
    slong nbz = (hasz ? 2 : 1);
    acb_t q4, q, u, v, w, t3, t1;
    slong k;
    int res, w_is_unit;

    acb_init(q4);
    acb_init(q);
    acb_init(u);
    acb_init(v);
    acb_init(w);
    acb_init(t3);
    acb_init(t1);

    for (k = 0; k < 2; k++)
    {
        prec = prec + FLINT_MAX(0, acb_theta_dist_addprec(&d[k]));
        prec = prec + FLINT_MAX(0, acb_theta_dist_addprec(&d0[k]));
    }

    /* compute q_{1/4}, q */
    acb_mul_2exp_si(q4, acb_mat_entry(tau, 0, 0), -2);
    acb_exp_pi_i(q4, q4, prec);
    acb_pow_ui(q, q4, 4, prec);

    /* compute v, w */
    acb_exp_pi_i(v, t, prec);
    acb_exp_pi_i(w, z, prec);
    w_is_unit = arb_is_zero(acb_imagref(z));

    acb_one(u);
    for (k = 0; k < nbt; k++)
    {
        if (k > 0)
        {
            acb_mul(u, u, v, prec);
        }
        acb_modular_theta_sum(t3, &th[2 * k + 1], &th[2 * k], t1,
            u, 1, q, 1, prec);
        acb_mul(&th[2 * k + 1], &th[2 * k + 1], q4, prec);
    }

    if (hasz)
    {
        acb_set(u, w);
        for (k = 0; k < nbt; k++)
        {
            if (k > 0)
            {
                acb_mul(u, u, v, prec);
            }
            acb_modular_theta_sum(t3, &th[2 * nbt + 2 * k + 1], &th[2 * nbt + 2 * k], t1,
                u, w_is_unit, q, 1, prec);
            acb_mul(&th[2 * nbt + 2 * k + 1], &th[2 * nbt + 2 * k + 1], q4, prec);
        }
    }
    res = _acb_vec_is_finite(th, 2 * nbz * nbt);

    acb_clear(q4);
    acb_clear(q);
    acb_clear(u);
    acb_clear(v);
    acb_clear(w);
    acb_clear(t3);
    acb_clear(t1);
    return res;
}

int acb_theta_ql_a0_naive(acb_ptr th, acb_srcptr t, acb_srcptr z, arb_srcptr d0,
    arb_srcptr d, const acb_mat_t tau, slong guard, slong prec)
{
    slong g = acb_mat_nrows(tau);
    if (g == 1)
    {
        return acb_theta_ql_a0_naive_g1(th, t, z, d0, d, tau, guard, prec);
    }
    else
    {
        return acb_theta_ql_a0_naive_gen(th, t, z, d0, d, tau, guard, prec);
    }
}
