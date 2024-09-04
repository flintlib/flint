/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_jet_00, state)
{
    slong iter;

    /* Test: agrees with jet_00_notransform */
    for (iter = 0; iter < 25 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 2);
        slong nb = 1 + n_randint(state, 3);
        slong ord = n_randint(state, 3);
        slong mprec = 100 + n_randint(state, 400);
        slong prec = mprec + 50;
        slong bits = n_randint(state, 4);
        slong nbth = acb_theta_jet_nb(ord, g);
        acb_mat_t tau;
        acb_ptr z;
        acb_ptr th, test;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(nb * g);
        th = _acb_vec_init(nb * nbth);
        test = _acb_vec_init(nb * nbth);

        /* Sample tau not too far from reduced domain */
        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_mat_scalar_mul_2exp_si(tau, tau, -1);
        for (k = 0; k < nb; k++)
        {
            acb_siegel_randtest_vec_reduced(z + k * g, state, tau, 0, prec);
        }
        _acb_vec_scalar_mul_2exp_si(z, z, nb * g, 1);

        flint_printf("g = %wd, prec = %wd, nb = %wd, ord = %wd, tau, z:\n", g, prec, nb, ord);
        acb_mat_printd(tau, 5);
        _acb_vec_printd(z, nb * g, 5);

        /* Call jet_00 at precision mprec, compare with jet_00_notransform */
        acb_theta_jet_00(th, z, nb, tau, ord, prec);
        acb_theta_jet_00_notransform(test, z, nb, tau, ord, prec);

        flint_printf("th, test:\n");
        _acb_vec_printd(th, nb * nbth, 5);
        _acb_vec_printd(test, nb * nbth, 5);

        if (!_acb_vec_overlaps(th, test, nb * nbth))
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, prec = %wd, nb = %wd, ord = %wd, tau, z:\n", g, prec, nb, ord);
            acb_mat_printd(tau, 5);
            _acb_vec_printd(z, nb * g, 5);
            flint_printf("th, test:\n");
            _acb_vec_printd(th, nb * nbth, 5);
            _acb_vec_printd(test, nb * nbth, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, nb * g);
        _acb_vec_clear(th, nb * nbth);
        _acb_vec_clear(test, nb * nbth);
    }

    TEST_FUNCTION_END(state);
}
