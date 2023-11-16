/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"

TEST_FUNCTION_START(acb_dirichlet_jacobi, state)
{
    slong prec = 128;
    ulong q;

    /* check Jacobi sums */
    for (q = 29 * 29; q > 1; q = q%2 ? 3*q+1 : q/2)
    {
        slong m1, m2;
        dirichlet_group_t G;
        dirichlet_char_t chi1, chi2;

        acb_t s1, s2;

        dirichlet_group_init(G, q);
        dirichlet_char_init(chi1, G);
        dirichlet_char_init(chi2, G);

        acb_init(s1);
        acb_init(s2);

        dirichlet_char_one(chi1, G);

        for (m1 = 0; m1 < 50; m1++)
        {

            dirichlet_char_one(chi2, G);

            for (m2 = 0; m2 < 50; m2++)
            {

                acb_dirichlet_jacobi_sum_naive(s1, G, chi1, chi2, prec);
                acb_dirichlet_jacobi_sum(s2, G, chi1, chi2, prec);

                if (!acb_overlaps(s1, s2))
                {
                    flint_printf("FAIL: J_%wu(%wu,%wu)",
                            q, chi1->n, chi2->n);
                    flint_printf("\nnaive ");
                    acb_printd(s1, 25);
                    flint_printf("\ndefault ");
                    acb_printd(s2, 25);
                    flint_printf("\n");
                    flint_printf("cond = %wu, %wu, %wu\n",
                            dirichlet_conductor_char(G, chi1),
                            dirichlet_conductor_char(G, chi2),
                            dirichlet_conductor_ui(G, nmod_mul(chi1->n, chi2->n, G->mod))
                            );
                    flint_abort();
                }
                if (dirichlet_char_next(chi2, G) < 0)
                    break;

            }

            if (dirichlet_char_next(chi1, G) < 0)
                break;

        }

        acb_clear(s1);
        acb_clear(s2);
        dirichlet_group_clear(G);
        dirichlet_char_clear(chi1);
        dirichlet_char_clear(chi2);
    }

    TEST_FUNCTION_END(state);
}
