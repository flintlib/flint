/*
    Copyright (C) 2016 Pascal Molin
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"

TEST_FUNCTION_START(acb_dirichlet_l_vec_hurwitz, state)
{
    slong iter;

    for (iter = 0; iter < 100 * 0.1 * flint_test_multiplier(); iter++)
    {
        ulong i, q;
        slong prec;
        dirichlet_group_t G;
        dirichlet_char_t chi;
        acb_t s, z;
        acb_ptr v;
        acb_dirichlet_hurwitz_precomp_t pre;

        prec = 50 + n_randint(state, 50);
        q = 1 + n_randint(state, 50);

        dirichlet_group_init(G, q);
        dirichlet_char_init(chi, G);

        acb_init(s);
        acb_one(s);
        acb_div_si(s, s, 2, prec);

        v = _acb_vec_init(G->phi_q);

        if (n_randint(state, 2))
        {
            ulong A, K, N;
            acb_dirichlet_hurwitz_precomp_choose_param(&A, &K, &N, s, G->phi_q, prec);
            acb_dirichlet_hurwitz_precomp_init(pre, s, acb_is_one(s), A, K, N, prec);
        }
        else
        {
            acb_dirichlet_hurwitz_precomp_init_num(pre, s, acb_is_one(s), G->phi_q, prec);
        }

        /* all at once */
        if (n_randint(state, 2))
            acb_dirichlet_l_vec_hurwitz(v, s, pre, G, prec);
        else
            acb_dirichlet_l_vec_hurwitz(v, s, NULL, G, prec);

        /* check with complete loop */
        i = 0;
        acb_init(z);
        dirichlet_char_one(chi, G);
        do {
            if (n_randint(state, 2))
                acb_dirichlet_l_hurwitz(z, s, pre, G, chi, prec);
            else
                acb_dirichlet_l_hurwitz(z, s, NULL, G, chi, prec);

            if (!acb_overlaps(z, v + i))
            {
                flint_printf("\n L value differ");
                flint_printf("\nL(1/2, %wu) single = ", chi->n);
                acb_printd(z, 20);
                flint_printf("\nL(1/2, %wu) multi = ", chi->n);
                acb_printd(v + i, 20);
                flint_printf("\n\n");
                acb_vec_printd(v, G->phi_q, 10);
                flint_printf("\n\n");
                flint_abort();
            }
            else if (acb_rel_accuracy_bits(z) < prec - 8
                        || acb_rel_accuracy_bits(v + i) < prec - 8)
            {
                    flint_printf("FAIL\n\n");
                    flint_printf("q = %wu\n", q);
                    flint_printf("\nL(1/2,chi_%wu(%wu,)) inaccurate\n", q, chi->n);
                    flint_printf("\nsingle =\n");
                    acb_printd(z, 30);
                    flint_printf("\ndft =\n");
                    acb_printd(v + i, 30);
                    flint_printf("\nerrors %wd & %wd [prec = %wd]\n",
                        acb_rel_accuracy_bits(z),
                        acb_rel_accuracy_bits(v + i), prec);
                    flint_abort();
             }

            i++;
        } while (dirichlet_char_next(chi, G) >= 0);

        acb_clear(s);
        _acb_vec_clear(v, G->phi_q);
        dirichlet_char_clear(chi);
        dirichlet_group_clear(G);

        acb_dirichlet_hurwitz_precomp_clear(pre);
    }

    TEST_FUNCTION_END(state);
}
