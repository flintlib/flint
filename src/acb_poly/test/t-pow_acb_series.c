/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"

TEST_FUNCTION_START(acb_poly_pow_acb_series, state)
{
    slong iter;

    /* compare with exp/log */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong prec, trunc;
        acb_poly_t f, g, h1, h2;
        acb_t c;

        prec = 2 + n_randint(state, 200);
        trunc = n_randint(state, 20);

        acb_poly_init(f);
        acb_poly_init(g);
        acb_poly_init(h1);
        acb_poly_init(h2);
        acb_init(c);

        /* generate binomials */
        if (n_randint(state, 20) == 0)
        {
            acb_randtest(c, state, prec, 10);
            acb_poly_set_coeff_acb(f, 0, c);
            acb_randtest(c, state, prec, 10);
            acb_poly_set_coeff_acb(f, 1 + n_randint(state, 20), c);
        }
        else
        {
            acb_poly_randtest(f, state, 1 + n_randint(state, 20), prec, 10);
        }

        acb_poly_randtest(h1, state, 1 + n_randint(state, 20), prec, 10);

        acb_randtest(c, state, prec, 10);
        acb_poly_set_acb(g, c);

        /* f^c */
        acb_poly_pow_acb_series(h1, f, c, trunc, prec);

        /* f^c = exp(c*log(f)) */
        acb_poly_log_series(h2, f, trunc, prec);
        acb_poly_mullow(h2, h2, g, trunc, prec);
        acb_poly_exp_series(h2, h2, trunc, prec);

        if (!acb_poly_overlaps(h1, h2))
        {
            flint_printf("FAIL\n\n");

            flint_printf("prec = %wd\n", prec);
            flint_printf("trunc = %wd\n", trunc);

            flint_printf("f = "); acb_poly_printd(f, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_printd(c, 15); flint_printf("\n\n");
            flint_printf("h1 = "); acb_poly_printd(h1, 15); flint_printf("\n\n");
            flint_printf("h2 = "); acb_poly_printd(h2, 15); flint_printf("\n\n");

            flint_abort();
        }

        acb_poly_pow_acb_series(f, f, c, trunc, prec);

        if (!acb_poly_overlaps(f, h1))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        acb_poly_clear(f);
        acb_poly_clear(g);
        acb_poly_clear(h1);
        acb_poly_clear(h2);
        acb_clear(c);
    }

    TEST_FUNCTION_END(state);
}
