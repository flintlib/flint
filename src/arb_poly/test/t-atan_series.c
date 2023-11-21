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
#include "arb_poly.h"

TEST_FUNCTION_START(arb_poly_atan_series, state)
{
    slong iter;

    for (iter = 0; iter < 100 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n, qbits, rbits1, rbits2;
        fmpq_poly_t A;
        arb_poly_t a, b, c, d;

        qbits = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 30);
        n = 1 + n_randint(state, 30);

        fmpq_poly_init(A);
        arb_poly_init(a);
        arb_poly_init(b);
        arb_poly_init(c);
        arb_poly_init(d);

        fmpq_poly_randtest(A, state, m, qbits);
        arb_poly_set_fmpq_poly(a, A, rbits1);

        arb_poly_randtest(b, state, 1 + n_randint(state, 30), rbits1, 5);

        arb_poly_atan_series(b, a, n, rbits2);

        /* Check 2 atan(x) = atan(2x/(1-x^2)) + C */
        arb_poly_mullow(c, a, a, n, rbits2);
        arb_poly_one(d);
        arb_poly_sub(c, d, c, rbits2);
        arb_poly_add(d, a, a, rbits2);

        if (arb_poly_length(c) != 0)
        {
            arb_poly_div_series(c, d, c, n, rbits2);
            arb_poly_atan_series(c, c, n, rbits2);
            arb_poly_add(d, b, b, rbits2);

            /* TODO: also check the first coefficient */
            arb_poly_set_coeff_si(c, 0, 0);
            arb_poly_set_coeff_si(d, 0, 0);

            if (!arb_poly_overlaps(c, d))
            {
                flint_printf("FAIL\n\n");
                flint_printf("bits2 = %wd\n", rbits2);

                flint_printf("A = "); fmpq_poly_print(A); flint_printf("\n\n");
                flint_printf("a = "); arb_poly_printd(a, 15); flint_printf("\n\n");
                flint_printf("b = "); arb_poly_printd(b, 15); flint_printf("\n\n");
                flint_printf("c = "); arb_poly_printd(c, 15); flint_printf("\n\n");
                flint_printf("d = "); arb_poly_printd(d, 15); flint_printf("\n\n");

                flint_abort();
            }
        }

        arb_poly_atan_series(a, a, n, rbits2);
        if (!arb_poly_equal(a, b))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        fmpq_poly_clear(A);
        arb_poly_clear(a);
        arb_poly_clear(b);
        arb_poly_clear(c);
        arb_poly_clear(d);
    }

    TEST_FUNCTION_END(state);
}
