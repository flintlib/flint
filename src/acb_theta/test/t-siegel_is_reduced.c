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

    flint_printf("siegel_is_reduced....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: correct values on some matrices */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 6);
        slong tol_exp = -10;
        slong j = n_randint(state, g);
        slong k = n_randint(state, g);
        acb_mat_t tau;

        acb_mat_init(tau, g, g);

        acb_mat_onei(tau);
        if (!acb_siegel_is_reduced(tau, tol_exp, prec))
        {
            flint_printf("FAIL\n");
            acb_mat_printd(tau, 5);
            flint_abort();
        }

        acb_add_si(acb_mat_entry(tau, j, k), acb_mat_entry(tau, j, k), 1);
        acb_set(acb_mat_entry(tau, j, k), acb_mat_entry(tau, k, j));
        if (acb_siegel_is_reduced(tau, tol_exp, prec))
        {
            flint_printf("FAIL\n");
            acb_mat_printd(tau, 5);
            flint_abort();
        }

        acb_mat_onei(tau);
        acb_mul_2exp_si(acb_mat_entry(tau, j, j), acb_mat_entry(tau, j, j), -1);
        if (acb_siegel_is_reduced(tau, tol_exp, prec))
        {
            flint_printf("FAIL\n");
            acb_mat_printd(tau, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
