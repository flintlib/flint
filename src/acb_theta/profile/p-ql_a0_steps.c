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
#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"

static int usage(char * argv[])
{
    flint_printf("usage: %s g pstep pmax\n", argv[0]);
    return 1;
}

int main(int argc, char * argv[])
{
    slong iter = 0;
    flint_rand_t state;
    slong g, n, pstep, pmax, prec;

    if (argc < 4)
    {
        return usage(argv);
    }

    g = atol(argv[1]);
    n = 1 << g;
    pstep = atol(argv[2]);
    pmax = atol(argv[3]);

    flint_randinit(state);

    /* Profile with different number of steps on reduced input */
    for (prec = pstep; prec <= pmax; prec += pstep)
    {
        int hast = iter % 2;
        int hasz = (iter % 4) / 2;
        slong nbt = (hast ? 3 : 1);
        slong nbz = (hasz ? 2 : 1);
        slong guard = 2 * ACB_THETA_LOW_PREC;
        slong lp = ACB_THETA_LOW_PREC;
        acb_mat_t tau;
        arb_mat_t cho;
        acb_ptr z, t, r;
        arb_ptr dist, dist0;
        arb_t test;
        slong nb_steps, split;
        slong k;
        int res = 0;
        iter++;

        acb_mat_init(tau, g, g);
        arb_mat_init(cho, g, g);
        z = _acb_vec_init(g);
        t = _acb_vec_init(g);
        r = _acb_vec_init(nbz * nbt * n);
        dist = _arb_vec_init(n);
        dist0 = _arb_vec_init(n);
        arb_init(test);

        while (!res)
        {
            acb_siegel_randtest_reduced(tau, state, prec, 4);
            arb_sub_si(test, acb_imagref(acb_mat_entry(tau, g - 1, g - 1)), 3, prec);
            res = arb_is_negative(test);
        }
        acb_siegel_cho(cho, tau, lp);

        for (k = 0; k < g; k++)
        {
            if (hasz)
            {
                acb_urandom(&z[k], state, prec);
            }
            if (hast)
            {
                arb_urandom(acb_realref(&t[k]), state, prec);
            }
        }
        acb_theta_dist_a0(dist, z, tau, lp);
        acb_theta_dist_a0(dist0, t, tau, lp);

        split = 0;
        nb_steps = acb_theta_ql_a0_nb_steps(cho, 0, prec);

        flint_printf("(g = %wd, prec = %wd, hast = %wd, hasz = %wd) ideal nb_steps: %wd, tau:\n",
            g, prec, hast, hasz, nb_steps);
        acb_mat_printd(tau, 2);

        for (k = -FLINT_MIN(nb_steps, 2); k <= 2; k++)
        {
            flint_printf("nb_steps = %wd: ", nb_steps + k);
            TIMEIT_START;
            res = acb_theta_ql_a0_steps(r, t, z, dist0, dist, tau, nb_steps + k, split,
                guard, prec, &acb_theta_ql_a0_naive);
            TIMEIT_STOP;
            if (res)
            {
                flint_printf("result (expected prec loss %wd):\n",
                    (guard + g) * (nb_steps + k));
                acb_printd(&r[0], 5);
                flint_printf("\n");
            }
            else
            {
                flint_printf("FAIL\n");
            }
        }
        flint_printf("\n");

        acb_mat_clear(tau);
        arb_mat_clear(cho);
        _acb_vec_clear(z, g);
        _acb_vec_clear(t, g);
        _acb_vec_clear(r, nbz * nbt * n);
        _arb_vec_clear(dist, n);
        _arb_vec_clear(dist0, n);
        arb_clear(test);
    }

    flint_randclear(state);
    flint_cleanup();
    return 0;
}
