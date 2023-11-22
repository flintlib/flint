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

TEST_FUNCTION_START(acb_theta_g2_sextic_chi5, state)
{
    slong iter;

    /* Test: agrees with sextic and chi5 */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 2;
        slong n = 1 << (2 * g);
        slong prec = 100 + n_randint(state, 100);
        slong bits = n_randint(state, 4);
        acb_mat_t tau;
        acb_ptr z, th;
        acb_poly_t s1, s2;
        acb_t c1, c2;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        th = _acb_vec_init(n);
        acb_poly_init(s1);
        acb_poly_init(s2);
        acb_init(c1);
        acb_init(c2);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_mat_scalar_mul_2exp_si(tau, tau, -2);

        acb_theta_g2_sextic_chi5(s1, c1, tau, prec);
        acb_theta_g2_sextic(s2, tau, prec);

        if (!acb_poly_overlaps(s1, s2))
        {
            flint_printf("FAIL (sextic)\n");
            acb_poly_printd(s1, 5);
            flint_printf("\n");
            acb_poly_printd(s2, 5);
            flint_printf("\n");
            flint_abort();
        }

        acb_theta_all(th, z, tau, 0, prec);
        acb_theta_g2_chi5(c2, th, prec);

        if (!acb_overlaps(c1, c2))
        {
            flint_printf("FAIL (chi5)\n");
            acb_printd(c1, 10);
            flint_printf("\n");
            acb_printd(c2, 10);
            flint_printf("\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(th, n);
        acb_poly_clear(s1);
        acb_poly_clear(s2);
        acb_clear(c1);
        acb_clear(c2);
    }

    TEST_FUNCTION_END(state);
}
