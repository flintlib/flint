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
#include "acb_poly.h"

TEST_FUNCTION_START(arb_poly_product_roots_complex, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong prec, len1, len2, i;
        arb_poly_t a;
        acb_poly_t b, c;
        arb_ptr r;
        acb_ptr r2, r3;

        prec = 2 + n_randint(state, 200);
        len1 = n_randint(state, 12);
        len2 = n_randint(state, 8);

        arb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);

        r = _arb_vec_init(len1);
        r2 = _acb_vec_init(len2);
        r3 = _acb_vec_init(len1 + 2 * len2);

        for (i = 0; i < len1; i++)
            arb_randtest(r + i, state, 1 + n_randint(state, 200), 3);
        for (i = 0; i < len2; i++)
            acb_randtest(r2 + i, state, 1 + n_randint(state, 200), 3);

        for (i = 0; i < len1; i++)
            acb_set_arb(r3 + i, r + i);
        for (i = 0; i < len2; i++)
            acb_set(r3 + len1 + i, r2 + i);
        for (i = 0; i < len2; i++)
            acb_conj(r3 + len1 + len2 + i, r2 + i);

        arb_poly_product_roots_complex(a, r, len1, r2, len2, prec);
        acb_poly_product_roots(b, r3, len1 + 2 * len2, prec);

        acb_poly_set_arb_poly(c, a);

        if (!acb_poly_overlaps(b, c))
        {
            flint_printf("FAIL\n\n");

            flint_printf("len1 = %wd, len2 = %wd\n\n", len1, len2);
            for (i = 0; i < len1; i++)
            {
                arb_printd(r + i, 15);
                flint_printf("\n");
            }
            for (i = 0; i < len2; i++)
            {
                acb_printd(r2 + i, 15);
                flint_printf("\n");
            }

            flint_printf("a = "); arb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_poly_printd(c, 15); flint_printf("\n\n");

            flint_abort();
        }

        arb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(c);
        _arb_vec_clear(r, len1);
        _acb_vec_clear(r2, len2);
        _acb_vec_clear(r3, len1 + 2 * len2);
    }

    TEST_FUNCTION_END(state);
}
