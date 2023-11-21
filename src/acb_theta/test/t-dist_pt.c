/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_dist_pt, state)
{
    slong iter;

    /* Test: symmetric using C * pt as offset */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 6);
        slong prec = ACB_THETA_LOW_PREC;
        slong bits = n_randint(state, 3);
        arb_mat_t C;
        arb_ptr v;
        arb_t d1, d2;
        slong * pt1;
        slong * pt2;
        slong k;

        arb_mat_init(C, g, g);
        v = _arb_vec_init(g);
        arb_init(d1);
        arb_init(d2);
        pt1 = flint_malloc(g * sizeof(slong));
        pt2 = flint_malloc(g * sizeof(slong));

        arb_mat_randtest_cho(C, state, prec, bits);
        arb_mat_transpose(C, C);

        for (k = 0; k < g; k++)
        {
            pt1[k] = n_randint(state, 100);
            pt2[k] = n_randint(state, 100);
        }

        for (k = 0; k < g; k++)
        {
            arb_set_si(&v[k], pt1[k]);
        }
        arb_mat_vector_mul_col(v, C, v, prec);
        acb_theta_dist_pt(d1, v, C, pt2, prec);

        for (k = 0; k < g; k++)
        {
            arb_set_si(&v[k], pt2[k]);
        }
        arb_mat_vector_mul_col(v, C, v, prec);
        acb_theta_dist_pt(d2, v, C, pt1, prec);

        if (!arb_overlaps(d1, d2))
        {
            flint_printf("FAIL\n");
            flint_printf("C:\n");
            arb_mat_printd(C, 5);
            flint_printf("distances:\n");
            arb_printd(d1, 10);
            flint_printf("\n");
            arb_printd(d2, 10);
            flint_printf("\n");
            flint_abort();
        }

        arb_mat_clear(C);
        _arb_vec_clear(v, g);
        arb_clear(d1);
        arb_clear(d2);
        flint_free(pt1);
        flint_free(pt2);
    }

    TEST_FUNCTION_END(state);
}
