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

    flint_printf("jet_all....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with jet_naive_all */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 2);
        slong ord = n_randint(state, 3);
        slong nb = acb_theta_jet_nb(ord, g);
        slong n2 = 1 << (2 * g);
        slong prec = 100 + n_randint(state, 400);
        acb_mat_t tau;
        acb_ptr z, dth, test;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        dth = _acb_vec_init(nb * n2);
        test = _acb_vec_init(nb * n2);

        /* Sample tau not too far from reduced domain */
        acb_siegel_randtest_nice(tau, state, prec);
        acb_mat_scalar_mul_2exp_si(tau, tau, -1);
        arb_urandom(acb_realref(acb_mat_entry(tau, 0, 0)), state, prec);
        for (k = 0; k < g; k++)
        {
            acb_urandom(z, state, prec);
        }

        acb_theta_jet_all(dth, z, tau, ord, prec);
        acb_theta_jet_naive_all(test, z, tau, ord, prec);

        if (!_acb_vec_overlaps(dth, test, nb * n2))
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, prec = %wd, ord = %wd, tau:\n", g, prec, ord);
            acb_mat_printd(tau, 5);
            flint_printf("th, test:\n");
            _acb_vec_printd(dth, nb * n2, 5);
            _acb_vec_printd(test, nb * n2, 5);
            flint_abort();
        }

        /* Sample garbage tau, computation should fail quickly */
        acb_mat_randtest(tau, state, prec, 10);
        acb_set_si(acb_mat_entry(tau, 0, 0), -1);
        acb_theta_jet_all(dth, z, tau, ord, prec);
        {
            if (acb_is_finite(&dth[0]))
            {
                flint_printf("FAIL (not infinite)\n");
                flint_abort();
            }
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(dth, nb * n2);
        _acb_vec_clear(test, nb * n2);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
