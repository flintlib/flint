/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_jet_all, state)
{
    slong iter;

    /* Test: agrees with jet_naive_all */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 2);
        slong ord = n_randint(state, 3);
        slong nb = acb_theta_jet_nb(ord, g);
        slong n2 = 1 << (2 * g);
        slong prec = 100 + n_randint(state, 400);
        slong bits = n_randint(state, 4);
        acb_mat_t tau;
        acb_ptr z, dth, test;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        dth = _acb_vec_init(nb * n2);
        test = _acb_vec_init(nb * n2);

        /* Sample tau not too far from reduced domain */
        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_mat_scalar_mul_2exp_si(tau, tau, -1);
        acb_siegel_randtest_vec(z, state, g, prec);

        /* Sometimes phony input too */
        if (n_randint(state, 20) == 0)
        {
            k = n_randint(state, g);
            arb_zero(acb_imagref(acb_mat_entry(tau, k, k)));
        }

        acb_theta_jet_all(dth, z, tau, ord, prec);
        acb_theta_jet_naive_all(test, z, tau, ord, prec);

        if (!_acb_vec_overlaps(dth, test, nb * n2))
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, prec = %wd, ord = %wd, tau, z:\n", g, prec, ord);
            _acb_vec_printd(z, g, 5);
            acb_mat_printd(tau, 5);
            flint_printf("th, test:\n");
            _acb_vec_printd(dth, nb * n2, 5);
            _acb_vec_printd(test, nb * n2, 5);
            flint_printf("difference:\n");
            _acb_vec_sub(test, dth, test, nb * n2, prec);
            _acb_vec_printd(test, nb * n2, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(dth, nb * n2);
        _acb_vec_clear(test, nb * n2);
    }

    TEST_FUNCTION_END(state);
}
