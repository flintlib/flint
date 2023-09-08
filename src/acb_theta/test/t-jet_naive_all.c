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

    flint_printf("jet_naive_all....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: first values match naive_all */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong prec = ACB_THETA_LOW_PREC + n_randint(state, 200);
        slong bits = n_randint(state, 4);
        slong ord = n_randint(state, 3);
        slong g = 1 + n_randint(state, 3);
        slong nb_z = n_randint(state, 3);
        slong n2 = 1 << (2 * g);
        slong nb = acb_theta_jet_nb(ord, g + 1);
        acb_mat_t tau;
        acb_ptr z, th, dth, test;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g * nb_z);
        th = _acb_vec_init(n2 * nb_z);
        dth = _acb_vec_init(nb * n2 * nb_z);
        test = _acb_vec_init(n2 * nb_z);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        for (k = 0; k < g * nb_z; k++)
        {
            acb_urandom(&z[k], state, prec);
        }

        acb_theta_jet_naive_all(dth, ord, z, nb_z, tau, prec);
        acb_theta_naive_all(th, z, nb_z, tau, prec);
        for (k = 0; k < n2 * nb_z; k++)
        {
            acb_set(&test[k], &dth[k * nb]);
        }

        if (!_acb_vec_overlaps(th, test, n2 * nb_z))
        {
            flint_printf("FAIL (overlap)\n");
            flint_printf("g = %wd, prec = %wd, ord = %wd, nb_z = %wd\n", g, prec, ord, nb_z);
            acb_mat_printd(tau, 5);
            _acb_vec_printd(z, g * nb_z, 5);
            flint_printf("naive_all:\n");
            _acb_vec_printd(th, n2 * nb_z, 5);
            flint_printf("test:\n");
            _acb_vec_printd(test, n2 * nb_z, 5);
            flint_printf("dth:\n");
            _acb_vec_printd(dth, n2 * nb * nb_z, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g * nb_z);
        _acb_vec_clear(th, n2 * nb_z);
        _acb_vec_clear(dth, nb * n2 * nb_z);
        _acb_vec_clear(test, n2 * nb_z);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
