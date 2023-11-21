/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_poly.h"
#include "acb_poly.h"

TEST_FUNCTION_START(acb_poly_atan_series, state)
{
    slong iter;

    for (iter = 0; iter < 100 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n, qbits, rbits1, rbits2;
        fmpq_poly_t A;
        acb_poly_t a, b, c, d;

        qbits = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 30);
        n = 1 + n_randint(state, 30);

        fmpq_poly_init(A);
        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);
        acb_poly_init(d);

        fmpq_poly_randtest(A, state, m, qbits);
        fmpq_poly_set_coeff_si(A, 0, 0);
        acb_poly_set_fmpq_poly(a, A, rbits1);

        acb_poly_atan_series(b, a, n, rbits2);

        /* Check 2 atan(x) = atan(2x/(1-x^2)) + C */
        acb_poly_mullow(c, a, a, n, rbits2);
        acb_poly_one(d);
        acb_poly_sub(c, d, c, rbits2);
        acb_poly_add(d, a, a, rbits2);

        if (acb_poly_length(c) != 0)
        {
            acb_poly_div_series(c, d, c, n, rbits2);
            acb_poly_atan_series(c, c, n, rbits2);
            acb_poly_add(d, b, b, rbits2);

            /* TODO: also check the first coefficient */
            acb_poly_set_coeff_si(c, 0, 0);
            acb_poly_set_coeff_si(d, 0, 0);

            if (!acb_poly_overlaps(c, d))
            {
                flint_printf("FAIL\n\n");
                flint_printf("bits2 = %wd\n", rbits2);

                flint_printf("A = "); fmpq_poly_print(A); flint_printf("\n\n");
                flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
                flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
                flint_printf("c = "); acb_poly_printd(c, 15); flint_printf("\n\n");
                flint_printf("d = "); acb_poly_printd(d, 15); flint_printf("\n\n");

                flint_abort();
            }
        }

        acb_poly_atan_series(a, a, n, rbits2);
        if (!acb_poly_equal(a, b))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        fmpq_poly_clear(A);
        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(c);
        acb_poly_clear(d);
    }

    TEST_FUNCTION_END(state);
}
