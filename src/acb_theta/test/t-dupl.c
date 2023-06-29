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

    flint_printf("dupl....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: compare with naive algorithm */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong prec = 100 + n_randint(state, 500);
        slong mag_bits = n_randint(state, 2);
        slong rad_exp = -5;
        slong n = 1 << g;

        acb_mat_t tau;
        acb_ptr z;
        arf_t rad;
        acb_ptr th;
        acb_ptr th_dupl;
        acb_ptr test;
        arb_t err;
        slong k;

        acb_mat_init(tau, g, g);
        arf_init(rad);
        z = _acb_vec_init(2 * g);
        th = _acb_vec_init(2 * n);
        th_dupl = _acb_vec_init(2 * n);
        test = _acb_vec_init(2 * n);
        arb_init(err);

        acb_siegel_randtest(tau, state, prec, mag_bits);
        arf_one(rad);
        arf_mul_2exp_si(rad, rad, rad_exp);
        for (k = 0; k < g; k++)
        {
            acb_urandom(&z[k], state, prec);
        }
        _acb_vec_scalar_mul_2exp_si(z, z, g, rad_exp);

        acb_theta_naive(th, z, 2, tau, prec);
        acb_mat_scalar_mul_2exp_si(tau, tau, 1);
        acb_theta_naive(th_dupl, z, 2, tau, prec);
        _acb_vec_sqr(th_dupl, th_dupl, 2 * n, prec);
        acb_theta_dupl(test, th, g, prec);

        /*
           flint_printf("g = %wd, prec = %wd, tau, z:\n", g, prec);
           acb_mat_printd(tau, 10);
           for (k = 0; k < g; k++)
           {
           acb_printd(&z[k], 10); flint_printf("\n");
           }
           flint_printf("theta:\n");
           for (k = 0; k < 2*n; k++)
           {
           acb_printd(&th[k], 10); flint_printf("\n");
           }
         */

        if (!_acb_vec_overlaps(test, th_dupl, 2 * n))
        {
            flint_printf("FAIL (overlap)\n");
            flint_printf("g = %wd, prec = %wd, tau, z:\n", g, prec);
            acb_mat_printd(tau, 10);
            _acb_vec_printd(z, g, 10);
            flint_printf("theta:\n");
            _acb_vec_printd(th, 2 * n, 10);
            flint_printf("dupl:\n");
            _acb_vec_printd(th_dupl, 2 * n, 10);
            flint_printf("test:\n");
            _acb_vec_printd(test, 2 * n, 10);
            fflush(stdout);
            flint_abort();
        }

        acb_mat_clear(tau);
        arf_clear(rad);
        _acb_vec_clear(z, 2 * g);
        _acb_vec_clear(th, 2 * n);
        _acb_vec_clear(th_dupl, 2 * n);
        _acb_vec_clear(test, 2 * n);
        arb_clear(err);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
