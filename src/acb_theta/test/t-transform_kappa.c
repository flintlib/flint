/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_transform_kappa, state)
{
    slong iter;

    /* Test: kappa and kappa2 agree */
    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong bits = n_randint(state, 4);
        slong prec = 200;
        fmpz_mat_t mat;
        fmpz_mat_t x;
        acb_mat_t tau;
        acb_t sqrtdet;
        slong kappa, kappa2;

        fmpz_mat_init(mat, 2 * g, 2 * g);
        fmpz_mat_init(x, 2, 2);
        acb_mat_init(tau, g, g);
        acb_init(sqrtdet);

        sp2gz_randtest(mat, state, bits);
        acb_siegel_randtest_reduced(tau, state, prec, bits);

        kappa = acb_theta_transform_kappa(sqrtdet, mat, tau, prec);
        kappa2 = acb_theta_transform_kappa2(mat);

        if (kappa % 4 != kappa2)
        {
            flint_printf("FAIL\n");
            flint_printf("tau, mat:\n");
            acb_mat_printd(tau, 5);
            fmpz_mat_print_pretty(mat);
            flint_printf("kappa = %wd, kappa2 = %wd, sqrtdet:\n", kappa, kappa2);
            acb_printd(sqrtdet, 5);
            flint_printf("\n");
            flint_abort();
        }

        fmpz_mat_clear(mat);
        fmpz_mat_clear(x);
        acb_mat_clear(tau);
        acb_clear(sqrtdet);
    }

    TEST_FUNCTION_END(state);
}

