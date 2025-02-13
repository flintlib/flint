/*
    Copyright (C) 2023 Jean Kieffer

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
    flint_printf("usage: %s g pmin pmax cst s exp\n", argv[0]);
    return 1;
}

int main(int argc, char * argv[])
{
    flint_rand_t state;
    timeit_t t0;
    slong g, n, pmin, pmax, prec, s, exp, reps;
    double tmin, tref, tcur;
    int cst;
    slong * pattern;
    slong * test_pattern;
    slong * delta;
    slong * best_delta;
    acb_mat_t tau;
    acb_ptr z, th;
    slong j, k, l, m;
    int useful;

    if (argc < 7)
    {
        return usage(argv);
    }

    g = atol(argv[1]);
    n = 1 << g;
    pmin = atol(argv[2]);
    pmax = atol(argv[3]);
    cst = atoi(argv[4]);
    s = atol(argv[5]);
    exp = atol(argv[6]);

    flint_rand_init(state);
    pattern = flint_malloc(g * sizeof(slong));
    test_pattern = flint_malloc(g * sizeof(slong));
    delta = flint_malloc(g * sizeof(slong));
    best_delta = flint_malloc(g * sizeof(slong));
    acb_mat_init(tau, g, g);
    z = _acb_vec_init(2 * g);
    th = _acb_vec_init(2 * n * n);

    prec = pmin;
    acb_siegel_randtest_compact(tau, state, 1, prec);
    for (j = g - s; j < g; j++)
    {
        for (k = 0; k < g; k++)
        {
            acb_mul_2exp_si(acb_mat_entry(tau, j, k), acb_mat_entry(tau, j, k), exp);
            acb_mul_2exp_si(acb_mat_entry(tau, k, j), acb_mat_entry(tau, k, j), exp);
        }
    }
    if (!cst)
    {
        acb_siegel_randtest_vec_reduced(z + g, state, 1, tau, 1, prec);
    }
    flint_printf("g = %wd, cst = %wd, s = %wd\n", g, cst, s);
    flint_printf("Values of tau, z:\n");
    acb_mat_printd(tau, 5);
    _acb_vec_printd(z + g, g, 5);

    for (prec = pmin; prec <= pmax; prec = ceil(1.5 * prec))
    {
        acb_theta_ql_nb_steps(pattern, tau, cst, prec);

        flint_printf("\nAt prec = %wd, suggested pattern:", prec);
        for (j = 0; j < g; j++)
        {
            flint_printf(" %wd", pattern[j]);
            best_delta[j] = 0;
        }
        flint_printf("\n");

        tref = 0;
        tmin = 0;
        m = 9 - 2 * g;
        for (k = 0; k < n_pow(m, g); k++)
        {
            l = k;
            for (j = 0; j < g; j++)
            {
                delta[j] = l % m;
                l = l / m;
                if (delta[j] > m/2)
                {
                    delta[j] -= m;
                }
            }

            useful = 1;
            for (j = 0; j < g; j++)
            {
                test_pattern[j] = pattern[j] + delta[j];
                if (test_pattern[j] < 0
                    || (j > 0 && test_pattern[j] > test_pattern[j - 1]))
                {
                    useful = 0;
                    break;
                }
            }
            if (!useful)
            {
                continue;
            }

            flint_printf("Testing pattern", delta);
            for (j = 0; j < g; j++)
            {
                flint_printf(" %wd", test_pattern[j]);
            }
            flint_printf(", ");
            if (k == 0)
            {
                /* Do it once prior to measuring ? */
                acb_theta_ql_exact(th, z, 2, tau, test_pattern, 1, 0, prec);
            }
            TIMEIT_REPEAT(t0, reps);
            acb_theta_ql_exact(th, z, 2, tau, test_pattern, 1, 0, prec);
            TIMEIT_END_REPEAT(t0, reps);

            tcur = ((double) t0->cpu) / reps;
            flint_printf("time: %f ms (%wd reps)\n", tcur, reps);

            if (k == 0)
            {
                tmin = tcur;
                tref = tcur;
            }
            else if (tcur < tmin)
            {
                tmin = tcur;
                for (j = 0; j < g; j++)
                {
                    best_delta[j] = delta[j];
                }
            }
        }

        if (tmin < tref)
        {
            flint_printf("\nAt prec = %wd, best pattern had t = %f ms:", prec, tmin);
            for (j = 0; j < g; j++)
            {
                flint_printf(" %wd", pattern[j] + best_delta[j]);
            }
            flint_printf("\nCompared to suggested pattern with t = %f ms\n", tref);
        }
        else
        {
            flint_printf("\nAt prec = %wd, suggested pattern was the best with t = %f ms\n",
                prec, tref);
        }
    }

    flint_rand_clear(state);
    flint_free(pattern);
    flint_free(test_pattern);
    flint_free(delta);
    flint_free(best_delta);
    acb_mat_clear(tau);
    _acb_vec_clear(z, g);
    _acb_vec_clear(th, 2 * n * n);

    flint_cleanup();
    return 0;
}
