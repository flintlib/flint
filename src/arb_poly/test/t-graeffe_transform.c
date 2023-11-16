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

TEST_FUNCTION_START(arb_poly_graeffe_transform, state)
{
    slong iter;
    arb_poly_t a, b, c;
    arb_ptr roots;
    arb_t leading;

    arb_poly_init(a);
    arb_poly_init(b);
    arb_poly_init(c);

    arb_init(leading);

    for (iter = 0; iter < 200 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong n, prec, i;

        n = n_randint(state, 20);
        prec = 2 + n_randint(state, 256);

        roots = _arb_vec_init(n);

        arb_randtest(leading, state, prec, n_randint(state, 16));

        for (i = 0; i < n; i++)
            arb_randtest(roots + i, state, prec, n_randint(state, 16));
        arb_poly_product_roots(a, roots, n, prec);
        arb_poly_scalar_mul(a, a, leading, prec);

        for (i = 0; i < n; i++)
            arb_sqr(roots + i, roots + i, prec);
        arb_sqr(leading, leading, prec);
        arb_poly_product_roots(c, roots, n, prec);
        arb_poly_scalar_mul(c, c, leading, prec);

        arb_poly_graeffe_transform(b, a, prec);
        if (!arb_poly_overlaps(b, c))
        {
            flint_printf("FAIL (overlap)\n\n");
            flint_printf("n = %wd, prec = %wd\n\n", n, prec);

            flint_printf("a: "); arb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b: "); arb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c: "); arb_poly_printd(c, 15); flint_printf("\n\n");

            flint_abort();
        }

        arb_poly_graeffe_transform(a, a, prec);
        if (!arb_poly_equal(a, b))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_printf("n = %wd, prec = %wd\n\n", n, prec);

            flint_printf("a: "); arb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b: "); arb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c: "); arb_poly_printd(c, 15); flint_printf("\n\n");

            flint_abort();
        }

        _arb_vec_clear(roots, n);
    }

    arb_poly_clear(a);
    arb_poly_clear(b);
    arb_poly_clear(c);

    arb_clear(leading);

    TEST_FUNCTION_END(state);
}
