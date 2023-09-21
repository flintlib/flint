/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"
#include "profiler.h"

static void
worker_dim1(acb_ptr th, acb_srcptr v1, acb_srcptr v2, const slong* precs, slong len,
    const acb_t cofactor, const slong* coords, slong ord, slong g, slong prec, slong fullprec)
{
    acb_t sum;

    acb_init(sum);

    acb_dot(sum, NULL, 0, v1, 1, v2, 1, len, prec);
    acb_addmul(th, sum, cofactor, fullprec);

    acb_clear(sum);
}

static void
acb_theta_naive_00_gen(acb_ptr th, acb_srcptr z, slong nb_z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_theta_eld_t E;
    acb_theta_precomp_t D;
    acb_ptr c;
    arb_ptr u;
    acb_ptr new_z;
    slong ord = 0;
    slong nb = 1;
    slong k;
    timeit_t txx;

    acb_theta_eld_init(E, g, g);
    acb_theta_precomp_init(D, nb_z, g);
    c = _acb_vec_init(nb_z);
    u = _arb_vec_init(nb_z);
    new_z = _acb_vec_init(g * nb_z);

    flint_printf("(naive_00) ellipsoid and precomp:\n");
    timeit_start(txx);
    acb_theta_eld_clear(E);
    acb_theta_precomp_clear(D);
    acb_theta_eld_init(E, g, g);
    acb_theta_precomp_init(D, nb_z, g);

    acb_theta_naive_ellipsoid(E, new_z, c, u, ord, z, nb_z, tau, prec);

    prec = acb_theta_naive_fullprec(E, prec);
    acb_theta_precomp_set(D, new_z, tau, E, prec);
    timeit_stop(txx);
      flint_printf("%wd ms\n", txx->cpu);
    

    flint_printf("(naive_00) worker:\n");
      timeit_start(txx);
    
    for (k = 0; k < nb_z; k++)
    {
        acb_theta_naive_worker(&th[k], nb, &c[k], &u[k], E, D, k, ord, prec, worker_dim1);
    }
    timeit_stop(txx);
      flint_printf("%wd ms\n", txx->cpu);

    acb_theta_eld_clear(E);
    acb_theta_precomp_clear(D);
    _acb_vec_clear(c, nb_z);
    _arb_vec_clear(u, nb_z);
    _acb_vec_clear(new_z, g * nb_z);
}

void
acb_theta_naive_00(acb_ptr th, acb_srcptr z, slong nb_z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong k;
    acb_ptr res;

    if (g == 1)
    {
        res = _acb_vec_init(4);
        for (k = 0; k < nb_z; k++)
        {
            acb_modular_theta(&res[0], &res[1], &res[2], &res[3], z + k * g,
                acb_mat_entry(tau, 0, 0), prec);
            acb_set(&th[k], &res[2]);
        }
        _acb_vec_clear(res, 4);
    }
    else
    {
        acb_theta_naive_00_gen(th, z, nb_z, tau, prec);
    }
}
