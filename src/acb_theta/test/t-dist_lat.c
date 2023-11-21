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
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_dist_lat, state)
{
    slong iter;

    /* Test: make ellipsoid to check it is indeed the minimal distance */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong prec = ACB_THETA_LOW_PREC;
        slong hprec = 200;
        slong bits = n_randint(state, 5);
        acb_mat_t tau;
        arb_mat_t C;
        acb_ptr z;
        arb_ptr v, y;
        arb_t d, test, x, s;
        arf_t R2;
        acb_theta_eld_t E;
        slong *pts;
        slong k;
        int r;

        acb_mat_init(tau, g, g);
        arb_mat_init(C, g, g);
        z = _acb_vec_init(g);
        v = _arb_vec_init(g);
        y = _arb_vec_init(g);
        arb_init(d);
        arb_init(test);
        arb_init(x);
        arb_init(s);
        acb_theta_eld_init(E, g, g);
        arf_init(R2);

        /* Get reduced C */
        acb_siegel_randtest_reduced(tau, state, hprec, bits);
        acb_siegel_randtest_vec(z, state, g, prec);

        _acb_vec_get_imag(v, z, g);
        acb_siegel_cho(C, tau, prec);

        acb_theta_dist_lat(d, v, C, prec);
        arb_get_ubound_arf(R2, d, prec);

        /* Test: ellipsoid has points and d is the minimum distance */
        r = acb_theta_eld_set(E, C, R2, v);

        if (r)
        {
            if (acb_theta_eld_nb_pts(E) == 0)
            {
                flint_printf("FAIL (no points)\n");
                flint_printf("g = %wd, C:\n", g);
                arb_mat_printd(C, 10);
                flint_printf("offset:\n");
                _arb_vec_printn(v, g, 10, 0);
                flint_printf("\n");
                flint_printf("Distance: ");
                arf_printd(R2, 10);
                flint_printf("\n");
                flint_abort();
            }

            pts = flint_malloc(acb_theta_eld_nb_pts(E) * sizeof(slong) * g);
            acb_theta_eld_points(pts, E);

            arb_pos_inf(test);
            for (k = 0; k < acb_theta_eld_nb_pts(E); k++)
            {
                acb_theta_dist_pt(x, v, C, pts + k * g, prec);
                arb_min(test, test, x, prec);
            }

            if (!arb_overlaps(d, test))
            {
                flint_printf("FAIL (wrong distance)\n");
                flint_printf("g = %wd, C:\n", g);
                arb_mat_printd(C, 10);
                flint_printf("offset:\n");
                _arb_vec_printn(v, g, 10, 0);
                flint_printf("\n");
                flint_printf("Distance: ");
                arf_printd(R2, 10);
                flint_printf("\n");
                flint_abort();
            }

            flint_free(pts);
        }

        acb_mat_clear(tau);
        arb_mat_clear(C);
        _acb_vec_clear(z, g);
        _arb_vec_clear(v, g);
        _arb_vec_clear(y, g);
        arb_clear(d);
        arb_clear(test);
        arb_clear(x);
        arb_clear(s);
        acb_theta_eld_clear(E);
        arf_clear(R2);
    }

    TEST_FUNCTION_END(state);
}
