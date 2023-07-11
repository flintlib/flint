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

    flint_printf("naive....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: other naive functions agree with naive_all */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong nb = n_pow(2, g);
        acb_mat_t tau;
        acb_ptr z;
        acb_ptr th;
        acb_ptr th_test;
        slong prec = 20 + n_randint(state, 100);
        slong mag_bits = n_randint(state, 2);
        ulong ab = n_randint(state, nb * nb);
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        th = _acb_vec_init(nb);
        th_test = _acb_vec_init(nb * nb);

        acb_siegel_randtest_reduced(tau, state, prec, mag_bits);
        for (k = 0; k < g; k++)
        {
            acb_urandom(&z[k], state, prec);
        }
        
        acb_theta_naive_all(th_test, z, 1, tau, prec);
        
        acb_theta_naive(th, z, 1, tau, prec);
        if (!_acb_vec_overlaps(th, th_test, nb))
        {
            flint_printf("FAIL (naive)\n");
            flint_printf("g = %wd, prec = %wd, tau:\n", g, prec);
            acb_mat_printd(tau, 10);
            flint_printf("z:\n");
            _acb_vec_printd(z, g, 10);
            flint_printf("th, th_test:\n");
            _acb_vec_printd(th, nb, 10);
            _acb_vec_printd(th_test, nb * nb, 10);
            fflush(stdout);
            flint_abort();
        }

        acb_theta_naive_ind(&th[0], ab, z, 1, tau, prec);
        if (!acb_overlaps(&th[0], &th_test[ab]))
        {
            flint_printf("FAIL (naive_ind)\n");
            flint_abort();
        }

        acb_theta_get_a0(th_test, th_test, g);
        acb_theta_naive_a0(th, z, 1, tau, prec);
        if (!_acb_vec_overlaps(th, th_test, nb))
        {            
            flint_printf("FAIL (naive_a0)\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(th, nb);
        _acb_vec_clear(th_test, nb * nb);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
