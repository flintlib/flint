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

TEST_FUNCTION_START(acb_poly_binomial_transform_basecase, state)
{
    slong iter;

    for (iter = 0; iter < 500 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_poly_t a, b, c, d;
        slong j, n, prec;

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);
        acb_poly_init(d);

        n = n_randint(state, 20);
        prec = 2 + n_randint(state, 200);

        acb_poly_randtest(a, state, n, prec, 10);
        acb_poly_randtest(b, state, n, prec, 10);
        acb_poly_randtest(c, state, n, prec, 10);

        /* check self-inversion property */
        acb_poly_binomial_transform_basecase(b, a, n, prec);
        acb_poly_binomial_transform_basecase(c, b, n, prec);

        acb_poly_set(d, a);
        acb_poly_truncate(d, n);

        if (!acb_poly_contains(c, d))
        {
            flint_printf("FAIL (containment)\n\n");
            flint_printf("n = %wd, prec = %wd\n\n", n, prec);

            flint_printf("a: "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b: "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c: "); acb_poly_printd(c, 15); flint_printf("\n\n");
            flint_printf("d: "); acb_poly_printd(d, 15); flint_printf("\n\n");

            flint_abort();
        }

        acb_poly_set(d, a);
        acb_poly_binomial_transform_basecase(d, d, n, prec);
        if (!acb_poly_equal(d, b))
        {
            flint_printf("FAIL (aliasing)\n\n");

            flint_printf("a: "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b: "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("d: "); acb_poly_printd(d, 15); flint_printf("\n\n");

            flint_abort();
        }

        /* compare with power series operations */
        acb_poly_zero(d);
        for (j = 1; j < n; j++)
            acb_poly_set_coeff_si(d, j, -1);
        acb_poly_compose_series(c, a, d, n, prec);
        for (j = 0; j < n; j++)
            acb_poly_set_coeff_si(d, j, 1);
        acb_poly_mullow(c, c, d, n, prec);

        if (!acb_poly_overlaps(b, c))
        {
            flint_printf("FAIL (power series)\n\n");

            flint_printf("a: "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b: "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c: "); acb_poly_printd(c, 15); flint_printf("\n\n");

            flint_abort();
        }

        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(c);
        acb_poly_clear(d);
    }

    TEST_FUNCTION_END(state);
}
