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

    flint_printf("jet_naive_00....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: values match jet_naive_all */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong prec = ACB_THETA_LOW_PREC + n_randint(state, 100);
        slong bits = n_randint(state, 4);
        slong ord = n_randint(state, 4);
        slong g = 1 + n_randint(state, 3);
        slong n2 = 1 << (2 * g);
        slong nb = acb_theta_jet_nb(ord, g);
        acb_mat_t tau;
        acb_ptr z, dth, test;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        dth = _acb_vec_init(nb);
        test = _acb_vec_init(nb * n2);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        for (k = 0; k < g; k++)
        {
            acb_urandom(&z[k], state, prec);
        }

        acb_theta_jet_naive_00(dth, z, tau, ord, prec);
        acb_theta_jet_naive_all(test, z, tau, ord, prec);

        if (!_acb_vec_overlaps(dth, test, nb))
        {
            flint_printf("FAIL (overlap)\n");
            flint_printf("g = %wd, prec = %wd, ord = %wd\n", g, prec, ord);
            acb_mat_printd(tau, 5);
            _acb_vec_printd(z, g, 5);
            flint_printf("jet_naive_00:\n");
            _acb_vec_printd(dth, nb, 5);
            flint_printf("jet_naive_all:\n");
            _acb_vec_printd(test, nb, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(dth, nb);
        _acb_vec_clear(test, nb * n2);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
