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

    flint_printf("uql....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with naive_a0 for g <= 4 */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 1; /* + n_randint(state, 4);*/
        slong n = 1 << g;
        slong prec = 1000 + n_randint(state, 1000);
        slong nb_z = 1; /* + n_randint(state, 10);*/
        acb_mat_t tau, entry;
        acb_ptr z, th, test;
        slong k;
        
        acb_mat_init(tau, g, g);
        acb_mat_init(entry, 1, 1);
        z = _acb_vec_init(nb_z * g);
        th = _acb_vec_init(n * nb_z);
        test = _acb_vec_init(n * nb_z);

        /* In general, use direct algorithm */
        acb_siegel_randtest_nice(tau, state, prec);
        for (k = 0; k < nb_z * g; k++)
        {
            acb_urandom(&z[k], state, prec);
        }
        acb_theta_uql(th, z, nb_z, tau, prec);
        acb_theta_naive_a0(test, z, nb_z, tau, prec);

        if (!_acb_vec_overlaps(th, test, n)
            || !acb_is_finite(&th[0]))
        {
            flint_printf("FAIL (generic)\n");
            flint_printf("g = %wd, prec = %wd\n", g, prec);
            acb_mat_printd(tau, 10);
            _acb_vec_printd(th, n * nb_z, 10);
            flint_printf("\n");
            _acb_vec_printd(test, n * nb_z, 10);
            flint_printf("\n");
            flint_abort();
        }

        /* Construct example with uql_roots: tau diagonal, z = (1+tau)/2 */
        acb_mat_zero(tau);
        for (k = 0; k < g; k++)
        {
            acb_siegel_randtest_nice(entry, state, prec);
            acb_set(acb_mat_entry(tau, k, k), acb_mat_entry(entry, 0, 0));
            acb_add_si(&z[k], acb_mat_entry(tau, k, k), 1, prec);
            acb_mul_2exp_si(&z[k], &z[k], -1);
        }
        acb_theta_uql(th, z, nb_z, tau, prec);
        acb_theta_naive_a0(test, z, nb_z, tau, prec);

        if (!_acb_vec_overlaps(th, test, n)
            || !acb_contains_zero(&test[n-1])
            || !acb_is_finite(&th[0]))
        {
            flint_printf("FAIL (special)\n");
            flint_printf("g = %wd, prec = %wd\n", g, prec);
            acb_mat_printd(tau, 10);
            _acb_vec_printd(th, n * nb_z, 10);
            flint_printf("\n");
            _acb_vec_printd(test, n * nb_z, 10);
            flint_printf("\n");
            flint_abort();
        }
        
        acb_mat_clear(tau);
        acb_mat_clear(entry);
        _acb_vec_clear(z, nb_z * g);        
        _acb_vec_clear(th, n * nb_z);
        _acb_vec_clear(test, n * nb_z);
    }
    
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}

