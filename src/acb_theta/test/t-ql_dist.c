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
    
    flint_printf("ql_dist....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: make ellipsoid to check it is indeed the minimal distance */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 6);
        slong prec = 100;
        slong bits = n_randint(state, 5);
        slong lowprec = 32;
        acb_mat_t tau;
        arb_mat_t cho;
        arb_ptr offset, y;
        arb_t dist, test, x, s;
        arf_t R2;
        acb_theta_eld_t E;
        slong *pts;
        slong k, j;

        acb_mat_init(tau, g, g);
        arb_mat_init(cho, g, g);
        offset = _arb_vec_init(g);
        y = _arb_vec_init(g);
        arb_init(dist);
        arb_init(test);
        arb_init(x);
        arb_init(s);
        acb_theta_eld_init(E, g, g);
        arf_init(R2);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_mat_get_imag(cho, tau);
        arb_mat_cho(cho, cho, prec);
        arb_mat_transpose(cho, cho);
        for (k = 0; k < g; k++)
        {
            arb_randtest_precise(&offset[k], state, prec, bits);
        }

        acb_theta_ql_dist(dist, offset, cho, lowprec);
        arb_sqr(x, dist, lowprec);
        arb_get_ubound_arf(R2, x, lowprec);

        /* Test: ellipsoid must have points; distance is indeed the minimum
           distance */
        arb_mat_scalar_mul_2exp_si(cho, cho, -1);
        acb_theta_eld_fill(E, cho, R2, offset, lowprec);

        if (acb_theta_eld_nb_pts(E) == 0)
        {
            flint_printf("FAIL (no points)\n");
            flint_abort();            
        }

        pts = flint_malloc(acb_theta_eld_nb_pts(E) * sizeof(slong));
        acb_theta_eld_points(pts, E);
        
        arb_pos_inf(test);
        for (k = 0; k < acb_theta_eld_nb_pts(E); k++)
        {
            for (j = 0; j < g; j++)
            {
                arb_set_si(&y[j], pts[k * g + j]);
            }
            arb_mat_vector_mul_col(y, cho, y, prec);
            _arb_vec_add(y, y, offset, g, prec);

            arb_zero(x);
            for (j = 0; j < g; j++)
            {
                arb_sqr(s, &y[j], prec);
                arb_add(x, x, s, prec);
            }
            arb_min(test, test, x, prec);
        }

        if (!arb_overlaps(dist, test))
        {
            flint_printf("FAIL (wrong distance)\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        arb_mat_clear(cho);
        _arb_vec_clear(y, g);
        arb_clear(dist);
        arb_clear(test);
        arb_clear(x);
        arb_clear(s);
        acb_theta_eld_clear(E);
        arf_clear(R2);
        flint_free(pts);
    }
    
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
