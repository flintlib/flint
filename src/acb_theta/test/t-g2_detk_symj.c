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

    flint_printf("g2_detk_symj....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: chain rule */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 2;
        acb_mat_t c1, c2, c3;
        acb_poly_t s, r, t;
        slong k = n_randint(state, 10);
        slong j = n_randint(state, 10);
        slong bits = 2;
        slong prec = 100 + n_randint(state, 200);

        acb_mat_init(c1, g, g);
        acb_mat_init(c2, g, g);
        acb_mat_init(c3, g, g);
        acb_poly_init(s);
        acb_poly_init(r);
        acb_poly_init(t);

        acb_mat_randtest(c1, state, prec, bits);
        acb_mat_randtest(c2, state, prec, bits);
        acb_mat_mul(c3, c1, c2, prec);
        acb_poly_randtest(s, state, j + 1, prec, bits);

        acb_theta_g2_detk_symj(r, c2, s, k, j, prec);
        acb_theta_g2_detk_symj(r, c1, r, k, j, prec);
        acb_theta_g2_detk_symj(t, c3, s, k, j, prec);

        if (!acb_poly_overlaps(t, r))
        {
            flint_printf("FAIL\n");
            acb_mat_printd(c1, 5);
            acb_mat_printd(c2, 5);
            flint_printf("source:\n");
            acb_poly_printd(s, 5);
            flint_printf("\nvalues:\n");
            acb_poly_printd(r, 5);
            flint_printf("\n");
            acb_poly_printd(t, 5);
            flint_printf("\n");
            flint_abort();
        }

        acb_mat_clear(c1);
        acb_mat_clear(c2);
        acb_mat_clear(c3);
        acb_poly_clear(s);
        acb_poly_clear(r);
        acb_poly_clear(t);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
