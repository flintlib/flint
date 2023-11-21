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

TEST_FUNCTION_START(acb_theta_eld_points, state)
{
    slong iter;

    /* Test: all ellipsoid points must be within the box
       Then, generate random points:
       - points inside ellipsoid must appear in all_pts
       - points outside ellipsoid must have norm greater than R2 */
    for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong prec = ACB_THETA_LOW_PREC;
        slong mag_bits = n_randint(state, 2);
        acb_theta_eld_t E;
        arb_mat_t C;
        arf_t R2;
        arb_ptr v;
        slong k, j;
        slong try;
        slong *all_pts;
        slong *pt;
        int res;
        arb_mat_t vec;
        arb_t sqr, sum;

        acb_theta_eld_init(E, g, g);
        arb_mat_init(C, g, g);
        arf_init(R2);
        v = _arb_vec_init(g);
        pt = flint_malloc(g * sizeof(slong));
        arb_mat_init(vec, g, 1);
        arb_init(sqr);
        arb_init(sum);

        arb_mat_randtest_cho(C, state, prec, mag_bits);
        arb_mat_transpose(C, C);
        arb_randtest_positive(sqr, state, prec, mag_bits);
        arf_set(R2, arb_midref(sqr));
        arf_mul_si(R2, R2, 1 + n_randint(state, 10), prec, ARF_RND_UP);

        for (k = 0; k < g; k++)
        {
            arb_randtest_precise(&v[k], state, prec, mag_bits);
        }

        res = acb_theta_eld_set(E, C, R2, v);
        if (!res)
        {
            flint_printf("FAIL (ellipsoid)\n");
            flint_abort();
        }

        all_pts = flint_malloc(acb_theta_eld_nb_pts(E) * g * sizeof(slong));
        acb_theta_eld_points(all_pts, E);

        for (k = 0; k < acb_theta_eld_nb_pts(E); k++)
        {
            for (j = 0; j < g; j++)
            {
                if (FLINT_ABS(all_pts[k * g + j]) > acb_theta_eld_box(E, j))
                {
                    flint_printf("FAIL: point outside box\n");
                    flint_printf("\n");
                    flint_abort();
                }
            }
        }

        for (try = 0; try < 100; try++)
        {
            for (k = 0; k < g; k++)
            {
                pt[k] = n_randint(state, acb_theta_eld_box(E, k) + 1);
            }
            if (acb_theta_eld_contains(E, pt))
            {
                for (k = 0; k < acb_theta_eld_nb_pts(E); k++)
                {
                    res = 1;
                    for (j = 0; j < g; j++)
                    {
                        if (all_pts[k * g + j] != pt[j])
                        {
                            res = 0;
                            break;
                        }
                    }
                    if (res == 1)
                    {
                        break;
                    }
                }
                if (!res)
                {
                    flint_printf("FAIL: point not listed:\n");
                    for (j = 0; j < g; j++)
                    {
                        flint_printf("%wd ", pt[j]);
                    }
                    flint_abort();
                }
            }

            if (!acb_theta_eld_contains(E, pt))
            {
                arb_mat_zero(vec);
                for (k = 0; k < g; k++)
                {
                    arb_set_si(arb_mat_entry(vec, k, 0), pt[k]);
                }

                arb_mat_mul(vec, C, vec, prec);
                arb_zero(sum);
                for (k = 0; k < g; k++)
                {
                    arb_add(arb_mat_entry(vec, k, 0),
                            arb_mat_entry(vec, k, 0), &v[k], prec);
                    arb_sqr(sqr, arb_mat_entry(vec, k, 0), prec);
                    arb_add(sum, sum, sqr, prec);
                }
                arb_sub_arf(sum, sum, R2, prec);
                if (arb_is_negative(sum))
                {
                    flint_printf("FAIL: small point not in ellipsoid\n");
                    for (j = 0; j < g; j++)
                    {
                        flint_printf("%wd ", pt[j]);
                    }
                    flint_printf("\nCholesky:\n");
                    arb_mat_printd(C, 10);
                    flint_printf("Norm of point: ");
                    arb_printd(sum, 10);
                    flint_printf("\nCoordinates:\n");
                    for (j = 0; j < g; j++)
                    {
                        arb_printd(arb_mat_entry(vec, j, 0), 10);
                        flint_printf("\n");
                    }
                    flint_printf("Upper bound: ");
                    arf_printd(R2, 10);
                    flint_printf("\ntotal nb of points = %wd\n", acb_theta_eld_nb_pts(E));
                    flint_printf("Offset:\n");
                    for (j = 0; j < g; j++)
                    {
                        arb_printd(&v[j], 10);
                        flint_printf("\n");
                    }
                    flint_printf("Points:\n");
                    for (k = 0; k < acb_theta_eld_nb_pts(E); k++)
                    {
                        for (j = 0; j < g; j++)
                        {
                            flint_printf("%wd ", all_pts[k * g + j]);
                        }

                        flint_printf("\n");
                    }
                    flint_abort();
                }
            }
        }

        acb_theta_eld_clear(E);
        arb_mat_clear(C);
        arf_clear(R2);
        _arb_vec_clear(v, g);
        flint_free(all_pts);
        flint_free(pt);
        arb_mat_clear(vec);
        arb_clear(sqr);
        arb_clear(sum);
    }

    TEST_FUNCTION_END(state);
}
