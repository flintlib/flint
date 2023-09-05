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

    flint_printf("naive_radius....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: tail is small */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 10);
        slong ord = n_randint(state, 10);
        slong prec = ACB_THETA_LOW_PREC;
        slong bits = n_randint(state, 5);
        slong exp = 10 + n_randint(state, 100);
        arb_mat_t Y;
        arf_t R2, eps, t;
        arb_t bound;

        arb_mat_init(Y, g, g);
        arf_init(R2);
        arf_init(eps);
        arf_init(t);
        arb_init(bound);

        arb_mat_randtest_cho(Y, state, prec, bits);
        arb_mat_transpose(Y, Y);

        acb_theta_naive_radius(R2, eps, Y, ord, exp);
        acb_theta_naive_tail(bound, R2, Y, ord, prec);
        arb_get_lbound_arf(t, bound, prec);

        if (arf_cmp(t, eps) > 0)
        {
            flint_printf("FAIL\n");
            arb_mat_printd(Y, 5);
            flint_printf("exp = %wd, ord = %wd, eps, R2:\n", exp, ord);
            arf_printd(eps, 10);
            flint_printf("\n");
            arf_printd(R2, 10);
            flint_printf("\nbound:\n");
            arb_printd(bound, 10);
            flint_printf("\n");
            flint_abort();
        }

        arb_mat_clear(Y);
        arf_clear(R2);
        arf_clear(eps);
        arf_clear(t);
        arb_clear(bound);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
