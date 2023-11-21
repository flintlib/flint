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

TEST_FUNCTION_START(acb_poly_graeffe_transform, state)
{
    slong iter;
    acb_poly_t a, b, c;
    acb_ptr roots;
    acb_t leading;

    acb_poly_init(a);
    acb_poly_init(b);
    acb_poly_init(c);

    acb_init(leading);

    for (iter = 0; iter < 200 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong n, prec, i;

        n = n_randint(state, 20);
        prec = 2 + n_randint(state, 256);

        roots = _acb_vec_init(n);

        acb_randtest(leading, state, prec, n_randint(state, 16));

        for (i = 0; i < n; i++)
            acb_randtest(roots + i, state, prec, n_randint(state, 16));
        acb_poly_product_roots(a, roots, n, prec);
        acb_poly_scalar_mul(a, a, leading, prec);

        for (i = 0; i < n; i++)
            acb_sqr(roots + i, roots + i, prec);
        acb_sqr(leading, leading, prec);
        acb_poly_product_roots(c, roots, n, prec);
        acb_poly_scalar_mul(c, c, leading, prec);

        acb_poly_graeffe_transform(b, a, prec);
        if (!acb_poly_overlaps(b, c))
        {
            flint_printf("FAIL (overlap)\n\n");
            flint_printf("n = %wd, prec = %wd\n\n", n, prec);

            flint_printf("a: "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b: "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c: "); acb_poly_printd(c, 15); flint_printf("\n\n");

            flint_abort();
        }

        acb_poly_graeffe_transform(a, a, prec);
        if (!acb_poly_equal(a, b))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_printf("n = %wd, prec = %wd\n\n", n, prec);

            flint_printf("a: "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b: "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c: "); acb_poly_printd(c, 15); flint_printf("\n\n");

            flint_abort();
        }

        _acb_vec_clear(roots, n);
    }

    acb_poly_clear(a);
    acb_poly_clear(b);
    acb_poly_clear(c);

    acb_clear(leading);

    TEST_FUNCTION_END(state);
}
