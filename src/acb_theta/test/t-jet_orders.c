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

    flint_printf("jet_orders....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: get the right index */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 6);
        slong ord = n_randint(state, 6);
        slong nb = acb_theta_jet_nb(ord, g);
        slong *orders;
        slong i = n_randint(state, nb);
        slong test;
        slong j, k;

        orders = flint_malloc(nb * g * sizeof(slong));

        acb_theta_jet_orders(orders, ord, g);
        test = acb_theta_jet_index(orders + i * g, g);

        if (test != i)
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, ord = %wd, nb = %wd\n", g, ord, nb);
            flint_printf("orders:\n");
            for (j = 0; j < nb; j++)
            {
                for (k = 0; k < g; k++)
                {
                    flint_printf("%wd ", orders[j * g + k]);
                }
                flint_printf("\n");
            }
            flint_printf("i = %wd, test = %wd\n", i, test);
            flint_abort();
        }

        flint_free(orders);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
