/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_jet_one_notransform, state)
{
    slong iter;

    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 2);
        slong n2 = 1 << (2 * g);
        slong nb = n_randint(state, 3);
        slong ord = n_randint(state, 3);
        slong nbth = acb_theta_jet_nb(ord, g);
        slong prec = 100 + n_randint(state, 400);
        slong bits = n_randint(state, 4);
        ulong ab = n_randint(state, n2);
        acb_mat_t tau;
        acb_ptr z;
        acb_ptr th, all, test;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(nb * g);
        th = _acb_vec_init(nb * nbth);
        test = _acb_vec_init(nb * nbth);
        all = _acb_vec_init(n2 * nb * nbth);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        for (k = 0; k < nb; k++)
        {
            acb_siegel_randtest_vec_reduced(z + k * g, state, tau, 0, prec);
        }

        acb_theta_jet_one_notransform(th, z, nb, tau, ord, ab, prec);

        acb_theta_jet_all_notransform(all, z, nb, tau, ord, prec);
        for (k = 0; k < nb; k++)
        {
            _acb_vec_set(test + k * nbth, all + k * n2 * nbth + ab * nbth, nbth);
        }

        /* flint_printf("g = %wd, prec = %wd, nb = %wd, ab = %wd, ord = %wd tau, z:\n", g, prec, nb, ab, ord);
        acb_mat_printd(tau, 5);
        _acb_vec_printd(z, nb * g, 5);
        flint_printf("th, test:\n");
        _acb_vec_printd(th, nb * nbth, 5);
        _acb_vec_printd(test, nb * nbth, 5); */

        if (!_acb_vec_overlaps(th, test, nb * nbth))
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, prec = %wd, nb = %wd, ab = %wd, ord = %wd, tau, z:\n", g, prec, nb, ab, ord);
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
        _acb_vec_clear(all, n2 * nb * nbth);
    }

    TEST_FUNCTION_END(state);
}
