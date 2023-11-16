/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"

TEST_FUNCTION_START(acb_poly_root_bound_fujiwara, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_poly_t a;
        acb_ptr roots;
        acb_t t;
        mag_t mag1, mag2;
        slong i, deg, prec;

        prec = 10 + n_randint(state, 400);
        deg = n_randint(state, 10);

        acb_init(t);
        acb_poly_init(a);
        mag_init(mag1);
        mag_init(mag2);
        roots = _acb_vec_init(deg);

        for (i = 0; i < deg; i++)
            acb_randtest(roots + i, state, prec, 1 + n_randint(state, 20));

        acb_poly_product_roots(a, roots, deg, prec);
        acb_randtest(t, state, prec, 1 + n_randint(state, 20));
        acb_poly_scalar_mul(a, a, t, prec);

        acb_poly_root_bound_fujiwara(mag1, a);

        for (i = 0; i < deg; i++)
        {
            acb_get_mag(mag2, roots + i);

            /* acb_get_mag gives an upper bound which due to rounding
               could be larger than mag1, so we pick a slightly
               smaller number */
            mag_mul_ui(mag2, mag2, 10000);
            mag_div_ui(mag2, mag2, 10001);

            if (mag_cmp(mag2, mag1) > 0)
            {
                flint_printf("FAIL\n");
                flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
                flint_printf("root = "); acb_printd(roots + i, 15); flint_printf("\n\n");
                flint_printf("mag1 = "); mag_printd(mag1, 10); flint_printf("\n\n");
                flint_printf("mag2 = "); mag_printd(mag2, 10); flint_printf("\n\n");
                flint_abort();
            }
        }

        _acb_vec_clear(roots, deg);
        acb_clear(t);
        acb_poly_clear(a);
        mag_clear(mag1);
        mag_clear(mag2);
    }

    TEST_FUNCTION_END(state);
}
