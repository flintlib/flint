/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "profiler.h"
#include "acb_mat.h"
#include "acb_theta.h"

static int usage(char * argv[])
{
    flint_printf("usage: %s g nb_steps hasz\n", argv[0]);
    return 1;
}

int main(int argc, char * argv[])
{
    slong prec;
    acb_mat_t tau;
    acb_ptr th, z, t;
    arb_ptr d0, d;
    slong g, n, nb_steps, nbz, k;
    slong guard = 16;
    int hasz;

    if (argc < 4)
    {
        return usage(argv);
    }

    g = atol(argv[1]);
    n = 1 << (2 * g);
    nb_steps = atol(argv[2]);
    hasz = (int) atol(argv[3]);
    nbz = (hasz ? 2 : 1);

    acb_mat_init(tau, g, g);
    th = _acb_vec_init(n * nbz);
    z = _acb_vec_init(g);
    t = _acb_vec_init(g);
    d0 = _arb_vec_init(1 << g);
    d = _arb_vec_init(1 << g);

    acb_mat_onei(tau);
    for (k = 0; k < g - 1; k++)
    {
        acb_onei(acb_mat_entry(tau, k, k + 1));
        acb_mul_2exp_si(acb_mat_entry(tau, k, k + 1), acb_mat_entry(tau, k, k + 1), -2);
        acb_set(acb_mat_entry(tau, k + 1, k), acb_mat_entry(tau, k, k + 1));
    }

    prec = 32;
    if (hasz)
    {
        acb_set_si(z, 2);
        acb_sqrt(z, z, prec);
    }
    acb_theta_dist_a0(d0, t, tau, prec);
    acb_theta_dist_a0(d, z, tau, prec);

    for (k = 0; k < nb_steps; k++)
    {
        if (hasz)
        {
            acb_set_si(z, 2);
            acb_sqrt(z, z, prec);
        }

        flint_printf("prec = %wd, acb_theta_naive_all:\n", prec);
        TIMEIT_START;
        acb_theta_naive_all(th, z, 1, tau, prec);
        TIMEIT_STOP;
        acb_printd(&th[0], 5);
        flint_printf("\n");

        flint_printf("prec = %wd, acb_theta_ql_a0:\n", prec);
        TIMEIT_START;
        acb_theta_ql_a0(th, t, z, d0, d, tau, guard, prec);
        TIMEIT_STOP;
        acb_printd(&th[hasz * n], 5);
        flint_printf("\n");

        flint_printf("prec = %wd, acb_theta_ql_all:\n", prec);
        TIMEIT_START;
        acb_theta_ql_all(th, z, tau, 0, prec);
        TIMEIT_STOP;
        acb_printd(&th[0], 5);
        flint_printf("\n");

        flint_printf("prec = %wd, acb_theta_all:\n", prec);
        TIMEIT_START;
        acb_theta_all(th, z, tau, 0, prec);
        TIMEIT_STOP;
        acb_printd(&th[0], 5);
        flint_printf("\n\n");

        prec *= 2;
    }

    acb_mat_clear(tau);
    _acb_vec_clear(th, n * nbz);
    _arb_vec_clear(d0, 1 << g);
    _arb_vec_clear(d, 1 << g);
    _acb_vec_clear(z, g);
    _acb_vec_clear(t, g);
}
