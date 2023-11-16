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

TEST_FUNCTION_START(acb_poly_mullow, state)
{
    slong iter;

    /* compare with fmpq_poly */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong qbits1, qbits2, rbits1, rbits2, rbits3, trunc;
        fmpq_poly_t A, B, C;
        acb_poly_t a, b, c, d;

        qbits1 = 2 + n_randint(state, 200);
        qbits2 = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        rbits3 = 2 + n_randint(state, 200);
        trunc = n_randint(state, 10);

        fmpq_poly_init(A);
        fmpq_poly_init(B);
        fmpq_poly_init(C);

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);
        acb_poly_init(d);

        fmpq_poly_randtest(A, state, 1 + n_randint(state, 10), qbits1);
        fmpq_poly_randtest(B, state, 1 + n_randint(state, 10), qbits2);

        fmpq_poly_mullow(C, A, B, trunc);

        acb_poly_set_fmpq_poly(a, A, rbits1);
        acb_poly_set_fmpq_poly(b, B, rbits2);

        acb_poly_mullow(c, a, b, trunc, rbits3);

        if (!acb_poly_contains_fmpq_poly(c, C))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits3 = %wd\n", rbits3);
            flint_printf("trunc = %wd\n", trunc);

            flint_printf("A = "); fmpq_poly_print(A); flint_printf("\n\n");
            flint_printf("B = "); fmpq_poly_print(B); flint_printf("\n\n");
            flint_printf("C = "); fmpq_poly_print(C); flint_printf("\n\n");

            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_poly_printd(c, 15); flint_printf("\n\n");

            flint_abort();
        }

        acb_poly_set(d, a);
        acb_poly_mullow(d, d, b, trunc, rbits3);
        if (!acb_poly_equal(d, c))
        {
            flint_printf("FAIL (aliasing 1)\n\n");
            flint_abort();
        }

        acb_poly_set(d, b);
        acb_poly_mullow(d, a, d, trunc, rbits3);
        if (!acb_poly_equal(d, c))
        {
            flint_printf("FAIL (aliasing 2)\n\n");
            flint_abort();
        }

        /* test squaring */
        acb_poly_set(b, a);
        acb_poly_mullow(c, a, b, trunc, rbits3);
        acb_poly_mullow(d, a, a, trunc, rbits3);
        if (!acb_poly_overlaps(c, d))  /* not guaranteed to be identical */
        {
            flint_printf("FAIL (squaring)\n\n");

            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_poly_printd(c, 15); flint_printf("\n\n");

            flint_abort();
        }

        acb_poly_mullow(a, a, a, trunc, rbits3);
        if (!acb_poly_equal(d, a))
        {
            flint_printf("FAIL (aliasing, squaring)\n\n");

            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("d = "); acb_poly_printd(d, 15); flint_printf("\n\n");

            flint_abort();
        }

        fmpq_poly_clear(A);
        fmpq_poly_clear(B);
        fmpq_poly_clear(C);

        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(c);
        acb_poly_clear(d);
    }

    /* compare with classical */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong bits, trunc;
        acb_poly_t a, b, ab, ab2;

        bits = 2 + n_randint(state, 200);
        trunc = n_randint(state, 10);

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(ab);
        acb_poly_init(ab2);

        acb_poly_randtest(a, state, 1 + n_randint(state, 10), bits, 5);
        acb_poly_randtest(b, state, 1 + n_randint(state, 10), bits, 5);

        acb_poly_mullow_classical(ab, a, b, trunc, bits);
        acb_poly_mullow(ab2, a, b, trunc, bits);

        if (!acb_poly_overlaps(ab, ab2))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits = %wd\n", bits);
            flint_printf("trunc = %wd\n", trunc);

            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("ab = "); acb_poly_printd(ab, 15); flint_printf("\n\n");
            flint_printf("ab2 = "); acb_poly_printd(ab2, 15); flint_printf("\n\n");

            flint_abort();
        }

        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(ab);
        acb_poly_clear(ab2);
    }

    TEST_FUNCTION_END(state);
}
