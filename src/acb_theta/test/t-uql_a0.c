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

    /* Test: agrees with naive_a0 */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << g;
        slong prec = 200;
        slong bits = n_randint(state, 5);
        slong nb_z = 1 + n_randint(state, 2);
        acb_mat_t tau;
        acb_ptr z, th, test;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g * nb_z);
        th = _acb_vec_init(n * nb_z);
        test = _acb_vec_init(n * nb_z);

        /* Possibly generate large imaginary parts, but z is small */
        acb_siegel_randtest_reduced(tau, state, prec, bits);
        for (k = 0; k < nb_z * g; k++)
        {
            acb_urandom(&z[k], state, prec);
        }
        acb_theta_uql_a0(th, z, nb_z, tau, prec);
        acb_theta_naive_a0(test, z, nb_z, tau, prec);

        if (!_acb_vec_overlaps(th, test, n * nb_z)
            || !acb_is_finite(&th[0]))
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, prec = %wd\n", g, prec);
            acb_mat_printd(tau, 10);
            _acb_vec_printd(th, n * nb_z, 10);
            flint_printf("\n");
            _acb_vec_printd(test, n * nb_z, 10);
            flint_printf("\n");
            flint_abort();
        }
        
        acb_mat_clear(tau);
        _acb_vec_clear(z, g * nb_z);
        _acb_vec_clear(th, n * nb_z);
        _acb_vec_clear(test, n * nb_z);
    }
    
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
