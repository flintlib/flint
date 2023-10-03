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

    flint_printf("siegel_reduce_real....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: real part <= 1 in large precision */
    for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong prec = 100 + n_randint(state, 500);
        slong mag_bits = 1 + n_randint(state, 5);

        acb_mat_t tau;
        fmpz_mat_t m;
        arb_t test;

        acb_mat_init(tau, g, g);
        fmpz_mat_init(m, 2 * g, 2 * g);
        arb_init(test);

        acb_siegel_randtest(tau, state, prec, mag_bits);
        acb_siegel_reduce_real(m, tau);
        acb_siegel_transform(tau, m, tau, prec);

        arb_abs(test, acb_realref(acb_mat_entry(tau, 0, 0)));
        arb_sub_si(test, test, 1, prec);
        if (arb_is_positive(test))
        {
            flint_printf("FAIL (not reduced)\n");
            acb_mat_printd(tau, 10);
            fmpz_mat_print_pretty(m);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        acb_mat_clear(tau);
        fmpz_mat_clear(m);
        arb_clear(test);

    }
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
