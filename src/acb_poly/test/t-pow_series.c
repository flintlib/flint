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

TEST_FUNCTION_START(acb_poly_pow_series, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong rbits1, rbits2, len;
        acb_poly_t a, b, c, d;

        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        len = n_randint(state, 25);

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);
        acb_poly_init(d);

        if (n_randint(state, 4) == 0)
            acb_poly_randtest(a, state, 1, rbits1, 25);
        else
            acb_poly_randtest(a, state, 1 + n_randint(state, 20), rbits1, 5);

        if (n_randint(state, 4) == 0)
            acb_poly_randtest(b, state, 1, rbits1, 25);
        else
            acb_poly_randtest(b, state, 1 + n_randint(state, 20), rbits1, 5);

        acb_poly_randtest(c, state, 1 + n_randint(state, 20), rbits1, 5);

        acb_poly_pow_series(c, a, b, len, rbits2);

        /* a^b = exp(b log a) */
        acb_poly_log_series(d, a, len, rbits2);
        acb_poly_mullow(d, d, b, len, rbits2);
        acb_poly_exp_series(d, d, len, rbits2);

        if (!acb_poly_overlaps(c, d))
        {
            flint_printf("FAIL (iter %wd)\n\n", iter);
            flint_printf("bits2 = %wd\n", rbits2);
            flint_printf("len = %wd\n", len);

            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_poly_printd(c, 15); flint_printf("\n\n");
            flint_printf("d = "); acb_poly_printd(d, 15); flint_printf("\n\n");

            flint_abort();
        }

        /* check aliasing */
        if (iter < 5000)
        {
            acb_poly_set(d, a);
            acb_poly_pow_series(d, d, b, len, rbits2);

            if (!acb_poly_overlaps(c, d))
            {
                flint_printf("FAIL (aliasing 1)\n");
                flint_abort();
            }

            acb_poly_set(d, b);
            acb_poly_pow_series(d, a, d, len, rbits2);

            if (!acb_poly_overlaps(c, d))
            {
                flint_printf("FAIL (aliasing 2)\n");
                flint_abort();
            }
        }

        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(c);
        acb_poly_clear(d);
    }

    TEST_FUNCTION_END(state);
}
