/*
    Copyright (C) 2025 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <math.h>
#include "profiler.h"
#include "ulong_extras.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"

static int usage(char * argv[])
{
    flint_printf("usage: %s g nb_steps all skew\n", argv[0]);
    return 1;
}

int main(int argc, char * argv[])
{
    flint_rand_t state;
    slong g, n, prec, nb_steps, guard, skew;
    int all, res;
    slong nb = 2;
    acb_mat_t tau;
    acb_ptr zs, t;
    acb_ptr rts, rts_all;
    arb_ptr distances;
    slong * easy_steps;
    slong j;

    if (argc < 5)
    {
        return usage(argv);
    }

    g = atol(argv[1]);
    n = 1 << g;
    prec = 100000;
    nb_steps = atol(argv[2]);
    all = atoi(argv[3]);
    skew = atol(argv[4]);

    flint_rand_init(state);
    acb_mat_init(tau, g, g);
    zs = _acb_vec_init(2 * g);
    t = _acb_vec_init(g);
    rts = _acb_vec_init(6 * n * nb_steps);
    rts_all = _acb_vec_init(6 * n * n);
    distances = _arb_vec_init(2 * n);
    easy_steps = flint_malloc(2 * sizeof(slong));

    acb_siegel_randtest_compact(tau, state, 1, prec);
    for (j = 0; j < g; j++)
    {
        acb_mul_si(acb_mat_entry(tau, j, g - 1), acb_mat_entry(tau, j, g - 1), skew, prec);
        acb_mul_si(acb_mat_entry(tau, g - 1, j), acb_mat_entry(tau, g - 1, j), skew, prec);
    }
    acb_siegel_randtest_vec_reduced(zs + g, state, 1, tau, 1, prec);
    acb_theta_eld_distances(distances, zs, 2, tau, 32);

    flint_printf("g = %wd, nb_steps = %wd, tau, z:\n", g, nb_steps);
    acb_mat_printd(tau, 5);
    _acb_vec_printd(zs, 2 * g, 5);
    flint_printf("distances:\n");
    _arb_vec_printd(distances, 2 * n, 5);

    TIMEIT_START;
    res = acb_theta_ql_setup(rts, rts_all, t, &guard, easy_steps,
        zs, 2, tau, distances, nb_steps, all, prec);
    TIMEIT_STOP;

    flint_printf("res = %wd, guard = %wd, easy_steps:", res, guard);
    for (j = 0; j < nb; j++)
    {
        flint_printf(" %wd", easy_steps[j]);
    }
    flint_printf("\n");

    flint_rand_clear(state);
    acb_mat_clear(tau);
    _acb_vec_clear(zs, 2 * g);
    _acb_vec_clear(t, g);
    _acb_vec_clear(rts, 6 * n * nb_steps);
    _acb_vec_clear(rts_all, 6 * n * n);
    _arb_vec_clear(distances, 2 * n);
    flint_free(easy_steps);
    flint_cleanup();
}
