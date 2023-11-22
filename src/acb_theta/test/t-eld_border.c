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

TEST_FUNCTION_START(acb_theta_eld_border, state)
{
    slong iter;

    /* Test: border points are not contained in the ellipsoid,
       nor any children */
    for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong prec = ACB_THETA_LOW_PREC;
        slong mag_bits = n_randint(state, 2);
        acb_theta_eld_t E;
        arb_mat_t C;
        arb_t x;
        arf_t R2;
        arb_ptr v;
        slong k, j;
        slong *all_pts;
        int r;

        acb_theta_eld_init(E, g, g);
        arb_mat_init(C, g, g);
        arf_init(R2);
        v = _arb_vec_init(g);
        arb_init(x);

        arb_mat_randtest_cho(C, state, prec, mag_bits);
        arb_mat_transpose(C, C);
        arb_randtest_positive(x, state, prec, mag_bits);
        arf_set(R2, arb_midref(x));
        arf_mul_si(R2, R2, 1 + n_randint(state, 10), prec, ARF_RND_UP);
        for (k = 0; k < g; k++)
        {
            arb_randtest_precise(&v[k], state, prec, mag_bits);
        }

        r = acb_theta_eld_set(E, C, R2, v);
        if (!r)
        {
            flint_printf("FAIL (ellipsoid)\n");
            flint_abort();
        }

        all_pts = flint_malloc(acb_theta_eld_nb_border(E) * g * sizeof(slong));
        acb_theta_eld_border(all_pts, E);

        for (k = 0; k < acb_theta_eld_nb_border(E); k++)
        {
            for (j = 0; j < g; j++)
            {
                if (acb_theta_eld_contains(E, all_pts + k * g))
                {
                    flint_printf("FAIL: point inside ellipsoid\n");
                    flint_abort();
                }
            }

            for (j = 0; j < acb_theta_eld_nr(E); j++)
            {
                if (acb_theta_eld_contains(acb_theta_eld_rchild(E, j), all_pts + k * g))
                {
                    flint_printf("FAIL: point inside right child %wd\n", j);
                    flint_abort();
                }
            }
            for (j = 0; j < acb_theta_eld_nl(E); j++)
            {
                if (acb_theta_eld_contains(acb_theta_eld_lchild(E, j), all_pts + k * g))
                {
                    flint_printf("FAIL: point inside left child %wd\n", j);
                    flint_abort();
                }
            }
        }

        acb_theta_eld_clear(E);
        arb_mat_clear(C);
        arf_clear(R2);
        _arb_vec_clear(v, g);
        arb_clear(x);
        flint_free(all_pts);
    }

    TEST_FUNCTION_END(state);
}
