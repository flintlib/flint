/*
    Copyright (C) 2025 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "profiler.h"
#include "acb.h"
#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"


static int usage(char * argv[])
{
    flint_printf("usage: %s g prec all_a all_b tilde skew\n", argv[0]);
    return 1;
}

int main(int argc, char * argv[])
{
    flint_rand_t state;
    slong g, n, prec, skew, nbth, j;
    int all_a, all_b, tilde;
    acb_mat_t tau;
    arb_mat_t im;
    acb_ptr z, th;
    arb_t det;
    acb_theta_ctx_tau_t ctx_tau;
    acb_theta_ctx_z_t ctx_z;
    arb_ptr distances;
    slong lp = ACB_THETA_LOW_PREC;

    if (argc < 7)
    {
        return usage(argv);
    }

    g = atol(argv[1]);
    n = 1 << g;
    prec = atol(argv[2]);
    all_a = atol(argv[3]);
    all_b = atoi(argv[4]);
    tilde = atoi(argv[5]);
    skew = atol(argv[6]);
    nbth = (all_a ? n : 1) * (all_b ? n : 1);

    flint_rand_init(state);
    acb_mat_init(tau, g, g);
    arb_mat_init(im, g, g);
    z = _acb_vec_init(g);
    acb_theta_ctx_tau_init(ctx_tau, all_a, g);
    acb_theta_ctx_z_init(ctx_z, g);
    distances = _arb_vec_init(n);
    th = _acb_vec_init(nbth);
    arb_init(det);

    flint_printf("Siegel reduction: ");
    TIMEIT_START;
    acb_siegel_randtest_compact(tau, state, 1, prec + lp);
    TIMEIT_STOP;

    for (j = 0; j < g; j++)
    {
        acb_mul_si(acb_mat_entry(tau, j, g - 1), acb_mat_entry(tau, j, g - 1), skew, prec + lp);
        acb_mul_si(acb_mat_entry(tau, g - 1, j), acb_mat_entry(tau, g - 1, j), skew, prec + lp);
    }
    acb_siegel_randtest_vec_reduced(z, state, 1, tau, 1, prec + lp);
    acb_mat_get_imag(im, tau);
    arb_mat_det(det, im, prec + lp);

    flint_printf("g = %wd, prec = %wd, tau, z:\n", g, prec);
    acb_mat_printd(tau, 5);
    _acb_vec_printd(z, g, 5);
    flint_printf("det(im tau): ");
    arb_printd(det, 5);
    flint_printf("\nComputing distances: ");
    TIMEIT_START;
    acb_theta_eld_distances(distances, z, 1, tau, lp);
    TIMEIT_STOP;

    flint_printf("Setting contexts: ");
    TIMEIT_START;
    acb_theta_ctx_tau_set(ctx_tau, tau, prec + lp);
    acb_theta_ctx_z_set(ctx_z, z, ctx_tau, prec + lp);
    TIMEIT_STOP;

    flint_printf("Evaluating theta: ");
    TIMEIT_START;
    acb_theta_sum(th, ctx_z, 1, ctx_tau, distances, all_a, all_b, tilde, prec);
    TIMEIT_STOP;

    acb_printd(&th[0], 5);
    flint_printf("\n");
    acb_printd(&th[nbth - 1], 5);
    flint_printf("\n");

    flint_rand_clear(state);
    acb_mat_clear(tau);
    arb_mat_clear(im);
    _acb_vec_clear(z, g);
    acb_theta_ctx_tau_clear(ctx_tau);
    acb_theta_ctx_z_clear(ctx_z);
    _arb_vec_clear(distances, n);
    _acb_vec_clear(th, nbth);
    arb_clear(det);
    flint_cleanup();
}
