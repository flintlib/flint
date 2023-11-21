/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_poly.h"

TEST_FUNCTION_START(arb_poly_root_bound_fujiwara, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_poly_t a;
        arb_ptr roots;
        arb_t t;
        mag_t mag1, mag2;
        slong i, deg, prec;

        prec = 10 + n_randint(state, 400);
        deg = n_randint(state, 10);

        arb_init(t);
        arb_poly_init(a);
        mag_init(mag1);
        mag_init(mag2);
        roots = _arb_vec_init(deg);

        for (i = 0; i < deg; i++)
            arb_randtest(roots + i, state, prec, 1 + n_randint(state, 20));

        arb_poly_product_roots(a, roots, deg, prec);
        arb_randtest(t, state, prec, 1 + n_randint(state, 20));
        arb_poly_scalar_mul(a, a, t, prec);

        arb_poly_root_bound_fujiwara(mag1, a);

        for (i = 0; i < deg; i++)
        {
            arb_get_mag(mag2, roots + i);

            /* arb_get_mag gives an upper bound which due to rounding
               could be larger than mag1, so we pick a slightly
               smaller number */
            mag_mul_ui(mag2, mag2, 10000);
            mag_div_ui(mag2, mag2, 10001);

            if (mag_cmp(mag2, mag1) > 0)
            {
                flint_printf("FAIL\n");
                flint_printf("a = "); arb_poly_printd(a, 15); flint_printf("\n\n");
                flint_printf("root = "); arb_printd(roots + i, 15); flint_printf("\n\n");
                flint_printf("mag1 = "); mag_printd(mag1, 10); flint_printf("\n\n");
                flint_printf("mag2 = "); mag_printd(mag2, 10); flint_printf("\n\n");
                flint_abort();
            }
        }

        _arb_vec_clear(roots, deg);
        arb_clear(t);
        arb_poly_clear(a);
        mag_clear(mag1);
        mag_clear(mag2);
    }

    TEST_FUNCTION_END(state);
}
