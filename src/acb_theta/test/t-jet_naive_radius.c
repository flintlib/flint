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

TEST_FUNCTION_START(acb_theta_jet_naive_radius, state)
{
    slong iter;

    /* Test: sum of terms on border of ellipsoid must be less than bound */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong mprec = 50 + n_randint(state, 100);
        slong prec = mprec + 50;
        slong bits = n_randint(state, 4);
        slong ord = n_randint(state, 4);
        slong nb = acb_theta_jet_nb(ord, g);
        acb_theta_eld_t E;
        acb_mat_t tau;
        arb_mat_t C;
        arf_t R2, eps;
        acb_ptr z, new_z;
        arb_ptr v, a;
        acb_t c, term;
        arb_t u, abs, sum;
        slong nb_pts;
        slong * pts;
        slong * tups;
        slong j, k;
        int res;

        acb_mat_init(tau, g, g);
        arb_mat_init(C, g, g);
        arf_init(R2);
        arf_init(eps);
        acb_theta_eld_init(E, g, g);
        z = _acb_vec_init(g);
        new_z = _acb_vec_init(g);
        v = _arb_vec_init(g);
        a = _arb_vec_init(g);
        acb_init(c);
        acb_init(term);
        arb_init(u);
        arb_init(abs);
        arb_init(sum);
        tups = flint_malloc(g * nb * sizeof(slong));

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        for (k = 0; k < g; k++)
        {
            acb_randtest_precise(&z[k], state, prec, bits);
        }
        acb_siegel_cho(C, tau, prec);
        acb_theta_naive_reduce(v, new_z, a, c, u, z, 1, tau, prec);

        acb_theta_jet_naive_radius(R2, eps, C, v, ord, mprec);
        arb_mul_arf(u, u, eps, prec);

        /* Test: sum of terms on the border of ellipsoid is less than u */
        res = acb_theta_eld_set(E, C, R2, v);
        if (!res)
        {
            flint_printf("FAIL (ellipsoid)\n");
            flint_abort();
        }

        nb_pts = acb_theta_eld_nb_border(E);
        pts = flint_malloc(g * nb_pts * sizeof(slong));
        acb_theta_eld_border(pts, E);
        acb_theta_jet_tuples(tups, ord, g);

        for (j = 0; j < nb; j++)
        {
            arb_zero(sum);
            for (k = 0; k < nb_pts; k++)
            {
                acb_theta_naive_term(term, new_z, tau, tups + j * g, pts + k * g, prec);
                acb_abs(abs, term, prec);
                arb_add(sum, sum, abs, prec);
            }

            acb_abs(abs, c, prec);
            arb_mul(sum, sum, abs, prec);
            arb_sub(abs, sum, u, prec);

            if (arb_is_positive(abs))
            {
                flint_printf("FAIL\n");
                flint_printf("sum, bound:\n");
                arb_printd(sum, 10);
                flint_printf("\n");
                arb_printd(u, 10);
                flint_printf("\ntau:\n");
                acb_mat_printd(tau, 5);
                flint_printf("z:\n");
                _acb_vec_printd(z, g, 10);
                acb_theta_eld_print(E);
                flint_abort();
            }
        }

        acb_mat_clear(tau);
        arb_mat_clear(C);
        arf_clear(R2);
        arf_clear(eps);
        acb_theta_eld_clear(E);
        _acb_vec_clear(z, g);
        _acb_vec_clear(new_z, g);
        _arb_vec_clear(v, g);
        _arb_vec_clear(a, g);
        acb_clear(c);
        acb_clear(term);
        arb_clear(u);
        arb_clear(abs);
        arb_clear(sum);
        flint_free(pts);
        flint_free(tups);
    }

    TEST_FUNCTION_END(state);
}
