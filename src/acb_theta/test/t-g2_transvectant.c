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
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_g2_transvectant, state)
{
    slong iter;

    /* Test: (f,f)_6 = -3*a3^2 + 8*a2*a4 - 20*a1*a5 + 120*a0*a6 */
    for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
        slong prec = 200;
        slong bits = 2;
        acb_poly_t f, g;
        acb_t c, test;

        acb_poly_init(f);
        acb_poly_init(g);
        acb_init(c);
        acb_init(test);

        acb_poly_randtest(f, state, 6, prec, bits);
        acb_poly_set_coeff_si(f, 6, 1);

        acb_theta_g2_transvectant(g, f, f, 6, 6, 6, prec);

        if (acb_poly_degree(g) > 0)
        {
            flint_printf("FAIL (degree)\n");
            acb_poly_printd(f, 5);
            flint_printf("\n");
            acb_poly_printd(g, 5);
            flint_printf("\n");
            flint_abort();
        }

        acb_mul(c, acb_poly_get_coeff_ptr(f, 0), acb_poly_get_coeff_ptr(f, 6), prec);
        acb_addmul_si(test, c, 120, prec);
        acb_mul(c, acb_poly_get_coeff_ptr(f, 1), acb_poly_get_coeff_ptr(f, 5), prec);
        acb_addmul_si(test, c, -20, prec);
        acb_mul(c, acb_poly_get_coeff_ptr(f, 2), acb_poly_get_coeff_ptr(f, 4), prec);
        acb_addmul_si(test, c, 8, prec);
        acb_sqr(c, acb_poly_get_coeff_ptr(f, 3), prec);
        acb_addmul_si(test, c, -3, prec);

        acb_poly_get_coeff_acb(c, g, 0);
        acb_mul_si(c, c, 60, prec);

        if (!acb_overlaps(test, c))
        {
            flint_printf("FAIL (value)\n");
            acb_printd(test, 5);
            flint_printf("\n");
            acb_printd(c, 5);
            flint_printf("\n");
            flint_abort();
        }

        acb_poly_clear(f);
        acb_poly_clear(g);
        acb_clear(c);
        acb_clear(test);
    }

    TEST_FUNCTION_END(state);
}

