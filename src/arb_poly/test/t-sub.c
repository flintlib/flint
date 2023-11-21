/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_poly.h"
#include "arb_poly.h"

TEST_FUNCTION_START(arb_poly_sub, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong qbits1, qbits2, rbits1, rbits2, rbits3;
        fmpq_poly_t A, B, C;
        arb_poly_t a, b, c;

        qbits1 = 2 + n_randint(state, 200);
        qbits2 = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        rbits3 = 2 + n_randint(state, 200);

        fmpq_poly_init(A);
        fmpq_poly_init(B);
        fmpq_poly_init(C);

        arb_poly_init(a);
        arb_poly_init(b);
        arb_poly_init(c);

        fmpq_poly_randtest(A, state, 1 + n_randint(state, 10), qbits1);
        fmpq_poly_randtest(B, state, 1 + n_randint(state, 10), qbits2);
        fmpq_poly_sub(C, A, B);

        arb_poly_set_fmpq_poly(a, A, rbits1);
        arb_poly_set_fmpq_poly(b, B, rbits2);
        arb_poly_sub(c, a, b, rbits3);

        if (!arb_poly_contains_fmpq_poly(c, C))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits3 = %wd\n", rbits3);

            flint_printf("A = "); fmpq_poly_print(A); flint_printf("\n\n");
            flint_printf("B = "); fmpq_poly_print(B); flint_printf("\n\n");
            flint_printf("C = "); fmpq_poly_print(C); flint_printf("\n\n");

            flint_printf("a = "); arb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); arb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); arb_poly_printd(c, 15); flint_printf("\n\n");

            flint_abort();
        }

        fmpq_poly_clear(A);
        fmpq_poly_clear(B);
        fmpq_poly_clear(C);

        arb_poly_clear(a);
        arb_poly_clear(b);
        arb_poly_clear(c);
    }

    TEST_FUNCTION_END(state);
}
