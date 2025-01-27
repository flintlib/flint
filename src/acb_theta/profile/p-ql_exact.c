/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "profiler.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"

static int usage(char * argv[])
{
    flint_printf("usage: %s g prec cst s\n", argv[0]);
    return 1;
}

int main(int argc, char * argv[])
{
    flint_rand_t state;
    slong g, n, prec, s;
    int cst;
    slong * pattern;
    slong * test_pattern;
    slong delta;
    acb_mat_t tau;
    acb_ptr z, th;
    slong j, k;

    if (argc < 5)
    {
        return usage(argv);
    }

    g = atol(argv[1]);
    n = 1 << g;
    prec = atol(argv[2]);
    cst = atoi(argv[3]);
    s = atol(argv[4]);

    flint_rand_init(state);
    pattern = flint_malloc(g * sizeof(slong));
    test_pattern = flint_malloc(g * sizeof(slong));
    acb_mat_init(tau, g, g);
    z = _acb_vec_init(2 * g);
    th = _acb_vec_init(2 * n * n);

    acb_siegel_randtest_compact(tau, state, 1, prec);
    for (j = g - s; j < g; j++)
    {
        for (k = 0; k < g; k++)
        {
            acb_mul_2exp_si(acb_mat_entry(tau, j, k), acb_mat_entry(tau, j, k), 3);
            acb_mul_2exp_si(acb_mat_entry(tau, k, j), acb_mat_entry(tau, k, j), 3);
        }
    }
    if (!cst)
    {
        acb_siegel_randtest_vec_reduced(z + g, state, 1, tau, 1, prec);
    }
    acb_theta_ql_nb_steps(pattern, tau, cst, prec);

    flint_printf("g = %wd, prec = %wd, cst = %wd, s = %wd\n", g, prec, cst, s);
    flint_printf("Values of tau, z:\n");
    acb_mat_printd(tau, 5);
    _acb_vec_printd(z + g, g, 5);
    flint_printf("Suggested pattern:");
    for (j = 0; j < g; j++)
    {
        flint_printf(" %wd", pattern[j]);
    }
    flint_printf("\n\n");

    for (delta = -4; delta <= 2; delta++)
    {
        for (j = 0; j < g; j++)
        {
            test_pattern[j] = pattern[j];
        }
        for (j = (s > 0 ? g - s : 0); j < g; j++)
        {
            test_pattern[j] = FLINT_MAX(0, pattern[j] + delta);
        }

        flint_printf("delta = %wd, testing pattern:", delta);
        for (j = 0; j < g; j++)
        {
            flint_printf(" %wd", test_pattern[j]);
        }
        flint_printf("\n");
        TIMEIT_START;
        acb_theta_ql_exact(th, z, 2, tau, test_pattern, 1, 0, prec);
        TIMEIT_STOP;
        flint_printf("th[0], th[n]: ");
        acb_printd(&th[0], 5);
        flint_printf("\n");
        acb_printd(&th[n], 5);
        flint_printf("\n\n");
    }

    flint_rand_clear(state);
    flint_free(pattern);
    flint_free(test_pattern);
    acb_mat_clear(tau);
    _acb_vec_clear(z, g);
    _acb_vec_clear(th, 2 * n * n);

    flint_cleanup();
    return 0;
}
