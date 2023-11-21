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
#include "acb_poly.h"

TEST_FUNCTION_START(acb_poly_divrem, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n, qbits1, qbits2, rbits1, rbits2, rbits3;
        fmpq_poly_t A, B, Q, R;
        acb_poly_t a, b, q, r;

        qbits1 = 2 + n_randint(state, 200);
        qbits2 = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        rbits3 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 20);
        n = 1 + n_randint(state, 20);

        fmpq_poly_init(A);
        fmpq_poly_init(B);
        fmpq_poly_init(Q);
        fmpq_poly_init(R);

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(q);
        acb_poly_init(r);

        fmpq_poly_randtest(A, state, m, qbits1);
        fmpq_poly_randtest_not_zero(B, state, n, qbits2);

        fmpq_poly_divrem(Q, R, A, B);

        acb_poly_set_fmpq_poly(a, A, rbits1);
        acb_poly_set_fmpq_poly(b, B, rbits2);

        acb_poly_divrem(q, r, a, b, rbits3);

        if (!acb_poly_contains_fmpq_poly(q, Q) ||
             !acb_poly_contains_fmpq_poly(r, R))
        {
            flint_printf("FAIL\n\n");

            flint_printf("A = "); fmpq_poly_print(A); flint_printf("\n\n");
            flint_printf("B = "); fmpq_poly_print(B); flint_printf("\n\n");
            flint_printf("Q = "); fmpq_poly_print(Q); flint_printf("\n\n");
            flint_printf("R = "); fmpq_poly_print(R); flint_printf("\n\n");

            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("q = "); acb_poly_printd(q, 15); flint_printf("\n\n");
            flint_printf("r = "); acb_poly_printd(r, 15); flint_printf("\n\n");

            flint_abort();
        }

        acb_poly_divrem(a, r, a, b, rbits3);
        if (!acb_poly_equal(a, q))
        {
            flint_printf("FAIL (aliasing q, a)\n\n");
        }
        acb_poly_set_fmpq_poly(a, A, rbits1);

        acb_poly_divrem(b, r, a, b, rbits3);
        if (!acb_poly_equal(b, q))
        {
            flint_printf("FAIL (aliasing q, b)\n\n");
            flint_abort();
        }
        acb_poly_set_fmpq_poly(b, B, rbits2);

        acb_poly_divrem(q, a, a, b, rbits3);
        if (!acb_poly_equal(a, r))
        {
            flint_printf("FAIL (aliasing r, a)\n\n");
            flint_abort();
        }
        acb_poly_set_fmpq_poly(a, A, rbits1);

        acb_poly_divrem(q, b, a, b, rbits3);
        if (!acb_poly_equal(b, r))
        {
            flint_printf("FAIL (aliasing r, b)\n\n");
            flint_abort();
        }

        fmpq_poly_clear(A);
        fmpq_poly_clear(B);
        fmpq_poly_clear(Q);
        fmpq_poly_clear(R);

        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(q);
        acb_poly_clear(r);
    }

    TEST_FUNCTION_END(state);
}
