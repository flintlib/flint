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
    flint_printf("usage: %s g prec cstep cmax\n", argv[0]);
    return 1;
}

int main(int argc, char * argv[])
{
    flint_rand_t state;
    slong g, n, prec, c, cstep, cmax;

    if (argc < 5)
    {
        return usage(argv);
    }

    g = atol(argv[1]);
    n = 1 << g;
    prec = atol(argv[2]);
    cstep = atol(argv[3]);
    cmax = atol(argv[4]);

    flint_randinit(state);

    /* Profile with different splittings on reduced input */
    for (c = cstep; c <= cmax; c += cstep)
    {
        slong guard = 2 * ACB_THETA_LOW_PREC;
        slong lp = ACB_THETA_LOW_PREC;
        acb_mat_t tau, tau1;
        arb_mat_t cho;
        acb_ptr t, r1, r2, r3;
        arb_ptr dist0;
        arb_t test;
        slong nb_steps_1, nb_steps_2, split;
        slong j, k;
        int res = 0;

        acb_mat_init(tau, g, g);
        acb_mat_init(tau1, g, g);
        arb_mat_init(cho, g, g);
        t = _acb_vec_init(g);
        r1 = _acb_vec_init(n);
        r2 = _acb_vec_init(n);
        r3 = _acb_vec_init(n);
        dist0 = _arb_vec_init(n);
        arb_init(test);

        while (!res)
        {
            acb_siegel_randtest_reduced(tau, state, prec, 4);
            arb_sub_si(test, acb_imagref(acb_mat_entry(tau, g - 1, g - 1)), 3, prec);
            res = arb_is_negative(test);
        }

        for (split = 1; split < g; split++)
        {
            acb_mat_set(tau1, tau);
            for (j = split; j < g; j++)
            {
                for (k = split; k < g; k++)
                {
                    acb_mul_si(acb_mat_entry(tau1, j, k), acb_mat_entry(tau1, j, k),
                        c, prec);
                }
            }

            flint_printf("g = %wd, prec = %wd, c = %wd, split = %wd, matrix:\n",
                g, prec, c, split);
            acb_mat_printd(tau1, 2);

            acb_theta_dist_a0(dist0, t, tau1, lp);
            acb_siegel_cho(cho, tau1, lp);
            nb_steps_1 = acb_theta_ql_a0_nb_steps(cho, split, prec);
            nb_steps_2 = acb_theta_ql_a0_nb_steps(cho, 0, prec);

            flint_printf("time for split (nb_steps = %wd):\n", nb_steps_1);
            TIMEIT_START;
            res = acb_theta_ql_a0_steps(r1, t, t, dist0, dist0, tau1, nb_steps_1,
                split, guard, prec, &acb_theta_ql_a0);
            TIMEIT_STOP;

            if (res)
            {

                flint_printf("time for non-split (nb_steps = %wd):\n", nb_steps_2);
                TIMEIT_START;
                res = acb_theta_ql_a0_steps(r2, t, t, dist0, dist0, tau1, nb_steps_2,
                    0, guard, prec, &acb_theta_ql_a0);
                TIMEIT_STOP;
            }

            if (res)
            {
                flint_printf("time for ql_a0:\n");
                TIMEIT_START;
                res = acb_theta_ql_a0(r3, t, t, dist0, dist0, tau1, guard, prec);
                TIMEIT_STOP;
            }

            if (res)
            {
                flint_printf("result for split (expected prec loss %wd):\n",
                    (guard + g) * nb_steps_1);
                acb_printd(&r1[0], 5);
                flint_printf("\nresult for non-split (expected prec loss %wd):\n",
                    (guard + g) * nb_steps_2);
                acb_printd(&r2[0], 5);
                flint_printf("\nresult for ql_a0:\n");
                acb_printd(&r3[0], 5);
                flint_printf("\n\n");
            }
            else
            {
                flint_printf("FAIL\n\n");
            }
        }

        acb_mat_clear(tau);
        acb_mat_clear(tau1);
        arb_mat_clear(cho);
        _acb_vec_clear(t, g);
        _acb_vec_clear(r1, n);
        _acb_vec_clear(r2, n);
        _acb_vec_clear(r3, n);
        _arb_vec_clear(dist0, n);
        arb_clear(test);
    }

    flint_randclear(state);
    flint_cleanup();
    return 0;
}
