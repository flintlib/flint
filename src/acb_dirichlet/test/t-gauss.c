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

TEST_FUNCTION_START(acb_dirichlet_gauss, state)
{
    slong prec = 128;
    ulong q;

    /* check Gauss sums */
    for (q = 3; q < 250; q ++)
    {
        dirichlet_group_t G;
        dirichlet_char_t chi;

        acb_t s1, s2, s3, s4;

        dirichlet_group_init(G, q);
        dirichlet_char_init(chi, G);

        acb_init(s1);
        acb_init(s2);
        acb_init(s3);
        acb_init(s4);
        dirichlet_char_one(chi, G);

        while (1) {

            acb_dirichlet_gauss_sum_naive(s1, G, chi, prec);
            acb_dirichlet_gauss_sum(s2, G, chi, prec);
            acb_dirichlet_gauss_sum_factor(s3, G, chi, prec);
            if (dirichlet_conductor_char(G, chi) == G->q)
                acb_dirichlet_gauss_sum_theta(s4, G, chi, prec);
            else
                acb_set(s4, s1);

            if (!acb_overlaps(s1, s2)
                    || !acb_overlaps(s1, s3)
                    || !acb_overlaps(s1, s4))
            {
                flint_printf("FAIL: G(chi_%wu(%wu))\n\n", q, chi->n);
                flint_printf("\nnaive ");
                acb_printd(s1, 25);
                flint_printf("\ndefault ");
                acb_printd(s2, 25);
                flint_printf("\nfactor ");
                acb_printd(s3, 25);
                flint_printf("\ntheta ");
                acb_printd(s4, 25);
                flint_abort();
            }

            if (dirichlet_char_next(chi, G) < 0)
                break;
        }
        acb_clear(s1);
        acb_clear(s2);
        acb_clear(s3);
        acb_clear(s4);

        dirichlet_group_clear(G);
        dirichlet_char_clear(chi);
    }

    TEST_FUNCTION_END(state);
}
