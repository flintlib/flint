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

TEST_FUNCTION_START(acb_theta_jet_tuples, state)
{
    slong iter;

    /* Test: get the right index */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 6);
        slong ord = n_randint(state, 6);
        slong nb = acb_theta_jet_nb(ord, g);
        slong * tups;
        slong i = n_randint(state, nb);
        slong test;
        slong j, k;

        tups = flint_malloc(nb * g * sizeof(slong));

        acb_theta_jet_tuples(tups, ord, g);
        test = acb_theta_jet_index(tups + i * g, g);

        if (test != i)
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, ord = %wd, nb = %wd\n", g, ord, nb);
            flint_printf("tups:\n");
            for (j = 0; j < nb; j++)
            {
                for (k = 0; k < g; k++)
                {
                    flint_printf("%wd ", tups[j * g + k]);
                }
                flint_printf("\n");
            }
            flint_printf("i = %wd, test = %wd\n", i, test);
            flint_abort();
        }

        flint_free(tups);
    }

    TEST_FUNCTION_END(state);
}
