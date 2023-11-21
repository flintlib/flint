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

TEST_FUNCTION_START(acb_theta_siegel_reduce, state)
{
    slong iter;

    /* Test: mat is symplectic and image passes acb_siegel_is_reduced */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong prec = 100 + n_randint(state, 200);
        slong mag_bits = n_randint(state, 5);
        slong tol_exp = -10;
        int fail = iter % 10;
        acb_mat_t tau;
        acb_mat_t w;
        fmpz_mat_t mat;

        acb_mat_init(tau, g, g);
        acb_mat_init(w, g, g);
        fmpz_mat_init(mat, 2 * g, 2 * g);

        if (fail)
        {
            mag_bits = 100;
        }
        acb_siegel_randtest(tau, state, prec, mag_bits);
        acb_siegel_reduce(mat, tau, prec);

        if (!sp2gz_is_correct(mat))
        {
            flint_printf("FAIL (symplectic)\n");
            fmpz_mat_print(mat);
            flint_abort();
        }

        acb_siegel_transform(w, mat, tau, prec);

        if (!fail && !acb_siegel_is_reduced(w, tol_exp, prec))
        {
            flint_printf("FAIL (not reduced)\n");
            acb_mat_printd(tau, 10);
            fmpz_mat_print_pretty(mat);
            flint_printf("\n");
            acb_mat_printd(w, 10);
            flint_abort();
        }

        acb_mat_clear(tau);
        acb_mat_clear(w);
        fmpz_mat_clear(mat);
    }

    TEST_FUNCTION_END(state);
}
