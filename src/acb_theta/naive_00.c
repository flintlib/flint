/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_modular.h"
#include "acb_theta.h"

static void
worker(acb_ptr th, acb_srcptr v1, acb_srcptr v2, const slong * precs, slong len,
    const acb_t cofactor, const slong * coords, slong ord, slong g, slong prec, slong fullprec)
{
    acb_t sum;

    acb_init(sum);

    acb_dot(sum, NULL, 0, v1, 1, v2, 1, len, prec);
    acb_addmul(th, sum, cofactor, fullprec);

    acb_clear(sum);
}

static void
acb_theta_naive_00_gen(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_theta_eld_t E;
    arb_mat_t C;
    arf_t R2, eps;
    acb_ptr cs;
    arb_ptr v, as, us;
    acb_ptr new_zs;
    slong k;
    int b;

    acb_theta_eld_init(E, g, g);
    arb_mat_init(C, g, g);
    arf_init(R2);
    arf_init(eps);
    cs = _acb_vec_init(nb);
    us = _arb_vec_init(nb);
    as = _arb_vec_init(g * nb);
    v = _arb_vec_init(g);
    new_zs = _acb_vec_init(g * nb);

    acb_siegel_cho(C, tau, prec);
    acb_theta_naive_radius(R2, eps, C, 0, prec);
    acb_theta_naive_reduce(v, new_zs, as, cs, us, zs, nb, tau, prec);
    b = acb_theta_eld_set(E, C, R2, v);

    if (b)
    {
        acb_theta_naive_worker(th, 1, new_zs, nb, tau, E, 0, prec, worker);

        for (k = 0; k < nb; k++)
        {
            acb_mul(&th[k], &th[k], &cs[k], prec);
            arb_mul_arf(&us[k], &us[k], eps, prec);
            acb_add_error_arb(&th[k], &us[k]);
        }
    }
    else
    {
        for (k = 0; k < nb; k++)
        {
            acb_indeterminate(&th[k]);
        }
    }

    acb_theta_eld_clear(E);
    arb_mat_clear(C);
    arf_clear(R2);
    arf_clear(eps);
    _acb_vec_clear(cs, nb);
    _arb_vec_clear(us, nb);
    _arb_vec_clear(as, g * nb);
    _arb_vec_clear(v, g);
    _acb_vec_clear(new_zs, g * nb);
}

static void
acb_theta_naive_00_g1(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec)
{
    acb_t q, w;
    acb_ptr res;
    int w_is_unit;
    slong k;

    acb_init(q);
    acb_init(w);
    res = _acb_vec_init(4);

    acb_exp_pi_i(q, acb_mat_entry(tau, 0, 0), prec);

    for (k = 0; k < nb; k++)
    {
        acb_exp_pi_i(w, &zs[k], prec);
        w_is_unit = arb_is_zero(acb_imagref(&zs[k]));
        acb_modular_theta_sum(&res[0], &res[1], &res[2], &res[3],
            w, w_is_unit, q, 1, prec);
        acb_set(&th[k], &res[2]);
    }

    acb_clear(q);
    acb_clear(w);
    _acb_vec_clear(res, 4);
}

void
acb_theta_naive_00(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);

    if (g == 1)
    {
        acb_theta_naive_00_g1(th, zs, nb, tau, prec);
    }
    else
    {
        acb_theta_naive_00_gen(th, zs, nb, tau, prec);
    }
}
