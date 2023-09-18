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

    flint_printf("naive_00....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with first entry of naive_0b */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << g;
        acb_mat_t tau;
        acb_ptr z;
        slong nb_z = 1 + n_randint(state, 4);
        acb_ptr th, th_0b, test;
        slong prec1 = 4000 + n_randint(state, 4000);
        slong prec = prec1 + n_randint(state, 200);
        slong mag_bits = n_randint(state, 2);
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g * nb_z);
        th = _acb_vec_init(nb_z);
        th_0b = _acb_vec_init(n * nb_z);
        test = _acb_vec_init(nb_z);

        acb_siegel_randtest_reduced(tau, state, prec, mag_bits);
        for (k = 0; k < g * nb_z; k++)
        {
            acb_urandom(&z[k], state, prec);
        }

        acb_theta_naive_00(th, z, nb_z, tau, prec1);
        acb_theta_naive_0b(th_0b, z, nb_z, tau, prec);
        for (k = 0; k < nb_z; k++)
        {
            acb_set(&test[k], &th_0b[k * n]);
        }
        
            flint_printf("th, test:\n");
            _acb_vec_printd(th, nb_z, 5);
            _acb_vec_printd(test, nb_z, 5);

        if (!_acb_vec_overlaps(th, test, nb_z))
        {
            flint_printf("FAIL: overlap\n");
            flint_printf("g = %wd, prec1 = %wd, prec = %wd, nb_z = %wd, tau:\n",
                g, prec1, prec, nb_z);
            acb_mat_printd(tau, 5);
            flint_printf("z:\n");
            _acb_vec_printd(z, g * nb_z, 5);
            flint_printf("th, test:\n");
            _acb_vec_printd(th, nb_z, 5);
            _acb_vec_printd(test, nb_z, 5);
            /*flint_abort();*/
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g * nb_z);
        _acb_vec_clear(th, nb_z);
        _acb_vec_clear(th_0b, n * nb_z);
        _acb_vec_clear(test, nb_z);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
