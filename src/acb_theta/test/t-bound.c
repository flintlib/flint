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

    flint_printf("bound....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: value of theta should be less than bound */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong prec = 100 + n_randint(state, 500);
        slong mag_bits = n_randint(state, 2);
        slong n = 1 << (2 * g);
        acb_mat_t tau;
        acb_ptr z;
        acb_t err;
        arb_t rad;
        arb_t bound;
        acb_ptr th;
        arb_t abs;
        slong j, k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        acb_init(err);
        arb_init(rad);
        arb_init(bound);
        th = _acb_vec_init(n);
        arb_init(abs);

        acb_siegel_randtest_reduced(tau, state, prec, mag_bits);
        for (k = 0; k < g; k++)
        {
            acb_urandom(&z[k], state, prec);
        }
        
        acb_theta_bound(arb_midref(rad), arb_midref(bound), z, tau, prec);

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
            for (k = 0; k < g; k++)
            {
                acb_urandom(err, state, prec);
                acb_mul_arb(err, err, rad, prec);
                acb_add(&z[k], &z[k], err, prec);
            }
        }
        acb_theta_naive_all(th, z, 1, tau, prec);

        _acb_vec_ninf(abs, th, n, prec);
        if (arb_gt(abs, bound))
        {
            flint_printf("FAIL: theta value is too large\n");
            flint_printf("g = %wd, prec = %wd, tau, z in disk:\n", g, prec);
            acb_mat_printd(tau, 10);
            _acb_vec_printd(z, g, 10);
            flint_printf("rad: ");
            arb_printd(rad, 10);
            flint_printf("\n");
            flint_printf("bound: ");
            arb_printd(bound, 10);
            flint_printf("\n");
            flint_printf("theta:\n");
            _acb_vec_printd(th, n, 10);
            fflush(stdout);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        acb_clear(err);
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
