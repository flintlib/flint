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

    flint_printf("siegel_randtest_nice....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: acb_siegel_reduce returns the identity matrix */
    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong prec = 100 + n_randint(state, 500);

        acb_mat_t tau;
        fmpz_mat_t mat;

        acb_mat_init(tau, g, g);
        fmpz_mat_init(mat, 2 * g, 2 * g);

        acb_siegel_randtest_nice(tau, state, prec);
        acb_siegel_reduce(mat, tau, prec);

        if (!fmpz_mat_is_one(mat))
        {
            flint_printf("FAIL (not reduced)\n");
            fmpz_mat_print(mat);
            acb_mat_printd(tau, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
        fmpz_mat_clear(tau);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
