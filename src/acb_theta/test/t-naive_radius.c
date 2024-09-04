/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_naive_radius, state)
{
    slong iter;

    /* Test: sum of terms on border of ellipsoid must be less than bound */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong mprec = 50 + n_randint(state, 100);
        slong prec = mprec + 50;
        slong bits = n_randint(state, 4);
        acb_theta_eld_t E;
        acb_mat_t tau;
        arb_mat_t C, Yinv;
        arf_t R2, eps;
        acb_ptr z;
        arb_ptr y, v;
        acb_t term;
        arb_t u, pi, abs, sum;
        slong nb_pts;
        slong * pts;
        slong k;
        int res;

        acb_mat_init(tau, g, g);
        arb_mat_init(C, g, g);
        arb_mat_init(Yinv, g, g);
        arf_init(R2);
        arf_init(eps);
        acb_theta_eld_init(E, g, g);
        z = _acb_vec_init(g);
        y = _arb_vec_init(g);
        v = _arb_vec_init(g);
        arb_init(u);
        arb_init(pi);
        acb_init(term);
        arb_init(abs);
        arb_init(sum);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_siegel_randtest_vec_reduced(z, state, tau, 0, prec);

        acb_siegel_cho(C, tau, prec);
        acb_siegel_yinv(Yinv, tau, prec);
        _acb_vec_get_imag(y, z, g);
        arb_mat_vector_mul_col(v, Yinv, y, prec);
        arb_dot(u, NULL, 0, v, 1, y, 1, g, prec);
        arb_const_pi(pi, prec);
        arb_mul(u, u, pi, prec);
        arb_exp(u, u, prec);
        arb_mat_vector_mul_col(v, C, v, prec);

        acb_theta_naive_radius(R2, eps, C, 0, mprec);
        arb_mul_arf(u, u, eps, prec);

        /* Test: sum of terms on the border of ellipsoid is less than u */
        res = acb_theta_eld_set(E, C, R2, v);
        if (!res)
        {
            flint_printf("FAIL (ellipsoid)\n");
            acb_mat_printd(tau, 5);
            arb_mat_printd(C, 5);
            _acb_vec_printd(z, g, 5);
            _arb_vec_printd(v, g, 5);
            flint_abort();
        }

        nb_pts = acb_theta_eld_nb_border(E);
        pts = flint_malloc(g * nb_pts * sizeof(slong));
        acb_theta_eld_border(pts, E);

        arb_zero(sum);
        for (k = 0; k < nb_pts; k++)
        {
            acb_theta_naive_term(term, z, tau, NULL, pts + k * g, prec);
            acb_abs(abs, term, prec);
            arb_add(sum, sum, abs, prec);
        }

        if (arb_gt(sum, u))
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

        acb_mat_clear(tau);
        arb_mat_clear(C);
        arb_mat_clear(Yinv);
        arf_clear(R2);
        arf_clear(eps);
        acb_theta_eld_clear(E);
        _acb_vec_clear(z, g);
        _arb_vec_clear(v, g);
        _arb_vec_clear(y, g);
        arb_clear(u);
        acb_clear(term);
        arb_clear(pi);
        arb_clear(abs);
        arb_clear(sum);
        flint_free(pts);
    }

    TEST_FUNCTION_END(state);
}
