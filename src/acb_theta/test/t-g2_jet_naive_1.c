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

TEST_FUNCTION_START(acb_theta_g2_jet_naive_1, state)
{
    slong iter;

    /* Test: agrees with usual jet_naive at the right indices */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 2;
        slong n = 1 << (2 * g);
        slong nb = acb_theta_jet_nb(1, g + 1);
        slong prec = 100 + n_randint(state, 3000);
        slong bits = n_randint(state, 4);
        acb_mat_t tau;
        acb_ptr z, dth, test;
        slong k;
        int res;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        dth = _acb_vec_init(n * nb);
        test = _acb_vec_init(n * nb);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        if (iter % 10 == 0)
        {
            acb_zero(acb_mat_entry(tau, 0, 0));
        }

        acb_theta_g2_jet_naive_1(dth, tau, prec);
        acb_theta_jet_naive_all(test, z, tau, 1, prec);

        for (k = 0; k < n; k++)
        {
            if (acb_theta_char_is_even(k, 2))
            {
                res = acb_overlaps(&dth[3 * k], &test[3 * k]);
            }
            else
            {
                res = _acb_vec_overlaps(&dth[3 * k + 1], &test[3 * k + 1], 2);
            }

            if (!res)
            {
                flint_printf("FAIL (k = %wd)\n", k);
                acb_mat_printd(tau, 5);
                flint_printf("values:\n");
                _acb_vec_printd(&dth[3 * k], 3, 5);
                _acb_vec_printd(&test[3 * k], 3, 5);
                flint_abort();
            }
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(dth, n * nb);
        _acb_vec_clear(test, n * nb);
    }

    TEST_FUNCTION_END(state);
}
