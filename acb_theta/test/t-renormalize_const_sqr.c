/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int
main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("renormalize_const_sqr....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with naive algorithm */
    for (iter = 0; iter < 50 * arb_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong nb = 1 << g;
        acb_mat_t tau;
        acb_ptr th2, th2_test;
        acb_t scal;
        slong prec = 100 + n_randint(state, 1000);

        int res;
        slong k;

        acb_mat_init(tau, g, g);
        th2 = _acb_vec_init(nb);
        th2_test = _acb_vec_init(nb);
        acb_init(scal);

        acb_siegel_randtest_fund(tau, state, prec);
        acb_theta_naive_const(th2, tau, prec);
        for (k = 0; k < nb; k++)
        {
            acb_sqr(&th2[k], &th2[k], prec);
        }
        _acb_vec_scalar_div(th2_test, th2, nb, &th2[0], prec);
        acb_theta_renormalize_const_sqr(scal, th2_test, tau, prec);
        _acb_vec_scalar_mul(th2_test, th2_test, nb, scal, prec);

        res = 1;
        for (k = 0; k < nb; k++)
        {
            if (!acb_overlaps(&th2_test[k], &th2[k]))
                res = 0;
        }
        if (!res)
        {
            flint_printf("FAIL (overlap)\n");
            flint_printf("th[k], th_test[k]:\n");
            for (k = 0; k < nb; k++)
            {
                acb_printd(&th2[k], 10);
                flint_printf("\n");
                acb_printd(&th2_test[k], 10);
                flint_printf("\n\n");
            }
            fflush(stdout);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(th2, nb);
        _acb_vec_clear(th2_test, nb);
        acb_clear(scal);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
