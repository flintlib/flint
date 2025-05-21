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
    flint_printf("usage: %s g prec ord all\n", argv[0]);
    return 1;
}

int main(int argc, char * argv[])
{
    flint_rand_t state;
    slong g, n, ord, prec, nbth, nbjet;
    int all;
    acb_mat_t tau;
    acb_ptr z, th;
    slong * pattern;
    acb_theta_ctx_tau_t ctx_tau;
    acb_theta_ctx_z_t ctx_z;
    slong j;

    if (argc < 5)
    {
        return usage(argv);
    }

    g = atol(argv[1]);
    prec = atol(argv[2]);
    ord = atol(argv[3]);
    all = atoi(argv[4]);

    n = 1 << g;
    nbjet = acb_theta_jet_nb(ord, g);
    nbth = (all ? n * n : n);

    flint_rand_init(state);
    acb_mat_init(tau, g, g);
    z = _acb_vec_init(g);
    th = _acb_vec_init(nbth * nbjet);
    pattern = flint_malloc(g * sizeof(slong));
    acb_theta_ctx_tau_init(ctx_tau, 0, g);
    acb_theta_ctx_z_init(ctx_z, g);

    acb_siegel_randtest_compact(tau, state, 0, prec);
    acb_siegel_randtest_vec_reduced(z, state, 1, tau, 0, prec);
    acb_theta_ql_nb_steps(pattern, tau, 0, prec);

    flint_printf("g = %wd, prec = %wd, pattern:", g, prec);
    for (j = 0; j < g; j++)
    {
        flint_printf(" %wd", pattern[j]);
    }

    flint_printf("\nql_jet_fd: ");
    TIMEIT_START;
    acb_theta_ql_jet_fd(th, z, 1, tau, ord, all, prec);
    TIMEIT_STOP;
    acb_printd(&th[nbth * nbjet - 1], 5);
    flint_printf("\nsum_jet: ");
    TIMEIT_START;
    acb_theta_ctx_tau_set(ctx_tau, tau, prec + ACB_THETA_LOW_PREC);
    acb_theta_ctx_z_set(ctx_z, z, ctx_tau, prec + ACB_THETA_LOW_PREC);
    acb_theta_sum_jet(th, ctx_z, 1, ctx_tau, ord, 1, all, prec);
    TIMEIT_STOP;
    acb_printd(&th[nbth * nbjet - 1], 5);
    flint_printf("\n");

    flint_rand_clear(state);
    acb_mat_clear(tau);
    _acb_vec_clear(z, g);
    _acb_vec_clear(th, nbth * nbjet);
    flint_free(pattern);
    acb_theta_ctx_tau_clear(ctx_tau);
    acb_theta_ctx_z_clear(ctx_z);
    flint_cleanup();
    return 0;
}
