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

    flint_printf("transform_sqrtdet_new....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: square of sqrtdet is det */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        acb_mat_t tau;
        acb_t r, t;
        slong prec = 100 + n_randint(state, 200);
        slong mag_bits = n_randint(state, 4);

        acb_mat_init(tau, g, g);
        acb_init(r);
        acb_init(t);

        acb_siegel_randtest(tau, state, prec, mag_bits);
        acb_theta_transform_sqrtdet_new(r, tau, prec);
        acb_sqr(r, r, prec);
        acb_mat_det(t, tau, prec);

        if (!acb_overlaps(r, t))
        {
            flint_printf("FAIL\n");
            acb_printd(r, 10);
            flint_printf("\n");
            acb_printd(t, 10);
            flint_printf("\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        acb_clear(r);
        acb_clear(t);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
