/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int
main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("bound_const....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: value of theta should be less than bound */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong prec = 100 + n_randint(state, 100);
        slong mag_bits = n_randint(state, 2);
        slong n = 1 << (2 * g);
        acb_mat_t tau;
        acb_t err;
        arb_t rad;
        arb_t bound;
        acb_ptr th;
        arb_t abs;
        slong j, k;

        acb_mat_init(tau, g, g);
        arb_init(rad);
        arb_init(bound);
        th = _acb_vec_init(n);
        arb_init(abs);

        acb_siegel_randtest(tau, state, prec, mag_bits);
        acb_theta_bound_const(arb_midref(rad), arb_midref(bound), tau, prec);

        /*
           flint_printf("g = %wd, prec = %wd, tau:\n", g, prec);
           acb_mat_printd(tau, 10);
           flint_printf("rad: "); arf_printd(rad, 10); flint_printf("\n");
           flint_printf("bound: "); arf_printd(bound, 10); flint_printf("\n");
         */

        if (!arb_is_positive(rad) || !arb_is_finite(bound))
        {
            flint_printf("Warning: not finite\n");
        }
        else
        {
            for (j = 0; j < g; j++)
            {
                for (k = j; k < g; k++)
                {
                    acb_urandom(err, state, prec);
                    acb_mul_arb(err, err, rad, prec);
                    acb_add(acb_mat_entry(tau, j, k), acb_mat_entry(tau, j, k),
                        err, prec);
                    acb_set(acb_mat_entry(tau, k, j), acb_mat_entry(tau, j, k));
                }
            }
        }
        acb_theta_naive_all_const(th, tau, prec);

        _acb_vec_ninf(abs, th, n, prec);
        if (arb_gt(abs, bound))
        {
            flint_printf("FAIL: theta value is too large\n");
            fflush(stdout);
            flint_abort();
        }

        acb_mat_clear(tau);
        arb_clear(rad);
        arb_clear(bound);
        _acb_vec_clear(th, n);
        arb_clear(abs);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
