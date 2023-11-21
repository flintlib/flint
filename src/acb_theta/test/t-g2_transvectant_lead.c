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

TEST_FUNCTION_START(acb_theta_g2_transvectant_lead, state)
{
    slong iter;

    /* Test: matches leading coefficient of g2_transvectant */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong prec = 200;
        slong bits = 2;
        acb_poly_t f, g, h;
        acb_t c, t;
        slong m = n_randint(state, 10);
        slong n = n_randint(state, 10);
        slong k = n_randint(state, FLINT_MIN(m, n) + 1);

        acb_poly_init(f);
        acb_poly_init(g);
        acb_poly_init(h);
        acb_init(c);
        acb_init(t);

        acb_poly_randtest(f, state, m, prec, bits);
        acb_poly_set_coeff_si(f, m, 1);
        acb_poly_randtest(g, state, n, prec, bits);
        acb_poly_set_coeff_si(g, n, 1);

        acb_theta_g2_transvectant(h, f, g, m, n, k, prec);
        acb_poly_get_coeff_acb(t, h, m + n - 2 * k);
        acb_theta_g2_transvectant_lead(c, f, g, m, n, k, prec);

        if (!acb_overlaps(c, t))
        {
            flint_printf("FAIL\n");
            flint_printf("m = %wd, n = %wd, k = %wd, f, g:\n", m, n, k);
            acb_poly_printd(f, 5);
            acb_poly_printd(g, 5);
            flint_printf("lead, test:\n");
            acb_printd(c, 5);
            flint_printf("\n");
            acb_printd(t, 5);
            flint_printf("\n");
            flint_abort();
        }

        acb_poly_clear(f);
        acb_poly_clear(g);
        acb_poly_clear(h);
        acb_clear(c);
        acb_clear(t);
    }

    TEST_FUNCTION_END(state);
}
