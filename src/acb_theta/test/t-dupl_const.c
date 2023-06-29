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

    flint_printf("dupl_const....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: coincide with theta duplication */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong prec = 100 + n_randint(state, 500);
        slong mag_bits = n_randint(state, 2);
        slong n = 1 << g;

        acb_mat_t tau;
        acb_ptr th;
        acb_ptr th_dupl;
        acb_ptr test;
        arb_t err;

        acb_mat_init(tau, g, g);
        th = _acb_vec_init(n);
        th_dupl = _acb_vec_init(n);
        test = _acb_vec_init(n);
        arb_init(err);

        acb_siegel_randtest(tau, state, prec, mag_bits);
        acb_theta_naive_const(th, tau, prec);
        acb_mat_scalar_mul_2exp_si(tau, tau, 1);
        acb_theta_naive_const(th_dupl, tau, prec);
        _acb_vec_sqr(th_dupl, th_dupl, n, prec);
        acb_theta_agm_step_sqrt(test, th, g, prec);

        /*
           flint_printf("g = %wd, prec = %wd, tau:\n", g, prec);
           acb_mat_printd(tau, 10);
           flint_printf("theta:\n");
           for (k = 0; k < n; k++)
           {
           acb_printd(&th[k], 10); flint_printf("\n");
           }
         */
        
        if (!_acb_vec_overlaps(test, th_dupl, n))
        {
            flint_printf("FAIL (overlap)\n");
            flint_printf("g = %wd, prec = %wd, tau:\n", g, prec);
            acb_mat_printd(tau, 10);
            flint_printf("theta:\n");
            _acb_vec_printd(th, n, 10);
            flint_printf("dupl:\n");
            _acb_vec_printd(th_dupl, n, 10);
            flint_printf("test:\n");
            _acb_vec_printd(test, n, 10);
            fflush(stdout);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(th, n);
        _acb_vec_clear(th_dupl, n);
        _acb_vec_clear(test, n);
        arb_clear(err);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
