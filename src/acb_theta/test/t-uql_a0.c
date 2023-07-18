/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("uql_a0....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with naive_ind */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 2;
        slong n = 1 << g;
        slong prec = 200 + n_randint(state, 200);
        slong bits = n_randint(state, 5);
        acb_mat_t tau;
        acb_ptr z, th, test;
        slong k;
        ulong a;
        
        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        th = _acb_vec_init(n);
        test = _acb_vec_init(n);

        /* Possibly generate large imaginary parts for tau and z */
        acb_siegel_randtest_reduced(tau, state, prec, bits);
        if (iter % 2 == 0)
        {
            acb_mul_2exp_si(acb_mat_entry(tau, g - 1, g - 1),
                acb_mat_entry(tau, g - 1, g - 1), 10);
        }
        for (k = 0; k < g; k++)
        {
            acb_urandom(&z[k], state, prec);
        }
        acb_theta_uql_a0(th, z, tau, prec);
        for (a = 0; a < n; a++)
        {
            acb_theta_naive_ind(&test[a], a << g, z, 1, tau, prec);
        }

        if (!_acb_vec_overlaps(th, test, n)
            || !acb_is_finite(&th[0]))
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, prec = %wd, tau:\n", g, prec);
            acb_mat_printd(tau, 10);
            flint_printf("z:\n");
            _acb_vec_printd(z, g, 10);
            flint_printf("\nvalues:\n");
            _acb_vec_printd(th, n, 10);
            flint_printf("\n");
            _acb_vec_printd(test, n, 10);
            flint_printf("\n");
            flint_abort();
        }
        
        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(th, n);
        _acb_vec_clear(test, n);
    }
    
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
