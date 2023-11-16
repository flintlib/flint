/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"

TEST_FUNCTION_START(acb_poly_product_roots, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong prec, len, m, i;
        acb_poly_t a, b, c, d;
        acb_ptr r;

        prec = 2 + n_randint(state, 200);
        len = n_randint(state, 12);
        if (len > 0)
            m = n_randint(state, len);
        else
            m = 0;

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);
        acb_poly_init(d);
        r = _acb_vec_init(len);

        for (i = 0; i < len; i++)
            acb_randtest(r + i, state, 1 + n_randint(state, 200), 10);

        acb_poly_product_roots(a, r, len, prec);
        acb_poly_product_roots(b, r, m, prec);
        acb_poly_product_roots(c, r + m, len - m, prec);
        acb_poly_mul(d, b, c, prec);

        if (!acb_poly_overlaps(a, d))
        {
            flint_printf("FAIL\n\n");

            flint_printf("len = %wd, m = %wd\n\n", len, m);
            for (i = 0; i < len; i++)
            {
                acb_printd(r + i, 15);
                flint_printf("\n");
            }

            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_poly_printd(c, 15); flint_printf("\n\n");
            flint_printf("d = "); acb_poly_printd(d, 15); flint_printf("\n\n");

            flint_abort();
        }

        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(c);
        acb_poly_clear(d);
        _acb_vec_clear(r, len);
    }

    TEST_FUNCTION_END(state);
}
