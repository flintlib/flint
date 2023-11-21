/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_poly.h"

TEST_FUNCTION_START(arb_poly_product_roots, state)
{
    slong iter;

    for (iter = 0; iter < 3000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong prec, len, m, i;
        arb_poly_t a, b, c, d;
        arb_ptr r;

        prec = 2 + n_randint(state, 200);
        len = n_randint(state, 12);
        if (len > 0)
            m = n_randint(state, len);
        else
            m = 0;

        arb_poly_init(a);
        arb_poly_init(b);
        arb_poly_init(c);
        arb_poly_init(d);
        r = _arb_vec_init(len);

        for (i = 0; i < len; i++)
            arb_randtest(r + i, state, 1 + n_randint(state, 200), 10);

        arb_poly_product_roots(a, r, len, prec);
        arb_poly_product_roots(b, r, m, prec);
        arb_poly_product_roots(c, r + m, len - m, prec);
        arb_poly_mul(d, b, c, prec);

        if (!arb_poly_overlaps(a, d))
        {
            flint_printf("FAIL\n\n");

            flint_printf("len = %wd, m = %wd\n\n", len, m);
            for (i = 0; i < len; i++)
            {
                arb_printd(r + i, 15);
                flint_printf("\n");
            }

            flint_printf("a = "); arb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); arb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); arb_poly_printd(c, 15); flint_printf("\n\n");
            flint_printf("d = "); arb_poly_printd(d, 15); flint_printf("\n\n");

            flint_abort();
        }

        arb_poly_clear(a);
        arb_poly_clear(b);
        arb_poly_clear(c);
        arb_poly_clear(d);
        _arb_vec_clear(r, len);
    }

    TEST_FUNCTION_END(state);
}
