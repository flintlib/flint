/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_poly.h"

TEST_FUNCTION_START(arb_poly_pow_arb_series, state)
{
    slong iter;

    /* compare with exp/log */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong prec, trunc;
        arb_poly_t f, g, h1, h2;
        arb_t c;

        prec = 2 + n_randint(state, 200);
        trunc = n_randint(state, 20);

        arb_poly_init(f);
        arb_poly_init(g);
        arb_poly_init(h1);
        arb_poly_init(h2);
        arb_init(c);

        /* generate binomials */
        if (n_randint(state, 20) == 0)
        {
            arb_randtest(c, state, prec, 10);
            arb_poly_set_coeff_arb(f, 0, c);
            arb_randtest(c, state, prec, 10);
            arb_poly_set_coeff_arb(f, 1 + n_randint(state, 20), c);
        }
        else
        {
            arb_poly_randtest(f, state, 1 + n_randint(state, 20), prec, 10);
        }

        arb_poly_randtest(h1, state, 1 + n_randint(state, 20), prec, 10);

        arb_randtest(c, state, prec, 10);
        arb_poly_set_arb(g, c);

        /* f^c */
        arb_poly_pow_arb_series(h1, f, c, trunc, prec);

        /* f^c = exp(c*log(f)) */
        arb_poly_log_series(h2, f, trunc, prec);
        arb_poly_mullow(h2, h2, g, trunc, prec);
        arb_poly_exp_series(h2, h2, trunc, prec);

        if (!arb_poly_overlaps(h1, h2))
        {
            flint_printf("FAIL\n\n");

            flint_printf("prec = %wd\n", prec);
            flint_printf("trunc = %wd\n", trunc);

            flint_printf("f = "); arb_poly_printd(f, 15); flint_printf("\n\n");
            flint_printf("c = "); arb_printd(c, 15); flint_printf("\n\n");
            flint_printf("h1 = "); arb_poly_printd(h1, 15); flint_printf("\n\n");
            flint_printf("h2 = "); arb_poly_printd(h2, 15); flint_printf("\n\n");

            flint_abort();
        }

        arb_poly_pow_arb_series(f, f, c, trunc, prec);

        if (!arb_poly_overlaps(f, h1))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        arb_poly_clear(f);
        arb_poly_clear(g);
        arb_poly_clear(h1);
        arb_poly_clear(h2);
        arb_clear(c);
    }

    TEST_FUNCTION_END(state);
}
