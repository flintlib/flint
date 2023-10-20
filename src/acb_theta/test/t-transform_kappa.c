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

    flint_printf("transform_kappa....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: matches combination of transform_kappa and transform_sqrtdet */
    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong bits = n_randint(state, 10);
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
        acb_siegel_randtest_nice(tau, state, prec);

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

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}

