/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_g2_sextic, state)
{
    slong iter;

    /* Test: discriminant of sextic is chi10 */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        slong g = 2;
        slong n = 1 << (2 * g);
        slong prec = 100 + n_randint(state, 100);
        slong bits = n_randint(state, 4);
        acb_mat_t tau;
        acb_ptr z, th2, roots;
        acb_poly_t f;
        acb_t d, t;
        slong nb, k, j;

        if (iter % 20 == 0)
        {
            prec += 10000; /* necessary for full test coverage */
        }

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        th2 = _acb_vec_init(n);
        roots = _acb_vec_init(6);
        acb_poly_init(f);
        acb_init(d);
        acb_init(t);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_mat_scalar_mul_2exp_si(tau, tau, -2);

        acb_theta_g2_sextic(f, tau, prec);
        nb = acb_poly_find_roots(roots, f, NULL, 0, prec);

        if (nb == 6)
        {
            acb_one(d);
            for (k = 0; k < 6; k++)
            {
                for (j = k + 1; j < 6; j++)
                {
                    acb_sub(t, &roots[k], &roots[j], prec);
                    acb_mul(d, d, t, prec);
                }
            }
            acb_sqr(d, d, prec);
            acb_poly_get_coeff_acb(t, f, 6);
            acb_pow_ui(t, t, 10, prec);
            acb_mul(d, d, t, prec);
            acb_mul_2exp_si(d, d, -12);

            acb_theta_all(th2, z, tau, 1, prec);
            acb_theta_g2_chi10(t, th2, prec);

            if (!acb_overlaps(d, t))
            {
                flint_printf("FAIL\n");
                flint_printf("roots, discr, chi10:\n");
                _acb_vec_printd(roots, 6, 5);
                acb_printd(d, 5);
                flint_printf("\n");
                acb_printd(t, 5);
                flint_printf("\n");
                flint_abort();
            }
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(th2, n);
        _acb_vec_clear(roots, 6);
        acb_poly_clear(f);
        acb_clear(d);
        acb_clear(t);
    }

    TEST_FUNCTION_END(state);
}
