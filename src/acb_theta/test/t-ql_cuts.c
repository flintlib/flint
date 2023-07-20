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
    
    flint_printf("ql_cuts....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: construct matrix that has such cuts */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 6);
        slong prec = 100;
        slong nb_cuts = n_randint(state, g);
        slong* cuts;
        slong nb_test;
        slong* test;
        arb_mat_t cho;
        arb_t c;
        slong i, j, k;
        int fail = 0;

        arb_mat_init(cho, g, g);
        arb_init(c);
        cuts = flint_malloc(nb_cuts * sizeof(slong));
        test = flint_malloc(g * sizeof(slong));

        /* Generate a 'random' strictly increasing sequence of cuts */
        for (k = 0; k < nb_cuts; k++)
        {
            cuts[k] = k + 1;
        }
        for (k = 1; k < g - nb_cuts; k++)
        {
            i = n_randint(state, nb_cuts + 1);
            for (j = i; j < nb_cuts; j++)
            {
                cuts[j]++;
            }
        }

        /* Make matrix */
        arb_one(c);
        i = 0;
        for (k = 0; k < g; k++)
        {
            arb_set(arb_mat_entry(cho, k, k), c);
            if ((i < nb_cuts) && (k + 1 == cuts[i]))
            {
                arb_mul_si(c, c, ACB_THETA_QL_CUT + 1, prec);
                i += 1;
            }
            else
            {
                arb_mul_si(c, c, ACB_THETA_QL_CUT - 1, prec);
            }
        }
        for (j = 0; j < g; j++)
        {
            for (k = j + 1; k < g; k++)
            {
                arb_urandom(arb_mat_entry(cho, j, k), state, prec);
            }
        }

        /* Compare */
        nb_test = acb_theta_ql_cuts(test, cho, prec);
        if (nb_test != nb_cuts)
        {
            fail = 1;
        }
        else
        {
            for (i = 0; i < nb_test; i++)
            {
                if (cuts[i] != test[i])
                {
                    fail = 1;
                    break;
                }
            }
        }

        if (fail)
        {
            flint_printf("FAIL\n");
            flint_printf("nb_cuts = %wd, cuts:\n", nb_cuts);
            for (k = 0; k < nb_cuts; k++)
            {
                flint_printf("%wd ", cuts[k]);
            }
            flint_printf("\n");
            flint_printf("nb_test = %wd, test:\n", nb_test);
            for (k = 0; k < nb_cuts; k++)
            {
                flint_printf("%wd ", test[k]);
            }
            flint_printf("\ncho:\n");
            arb_mat_printd(cho, 5);
            flint_abort();
        }

        arb_mat_clear(cho);
        arb_clear(c);
        flint_free(cuts);
        flint_free(test);
    }
    
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
