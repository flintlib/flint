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

FLINT_DLL extern slong arb_poly_newton_exp_cutoff;

TEST_FUNCTION_START(arb_poly_exp_series, state)
{
    slong iter;

    /* compare with fmpq_poly */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n, qbits, rbits1, rbits2;
        fmpq_poly_t A, B;
        arb_poly_t a, b;

        qbits = 2 + n_randint(state, 100);
        rbits1 = 2 + n_randint(state, 100);
        rbits2 = 2 + n_randint(state, 100);

        m = 1 + n_randint(state, 20);
        n = 1 + n_randint(state, 20);

        fmpq_poly_init(A);
        fmpq_poly_init(B);

        arb_poly_init(a);
        arb_poly_init(b);

        fmpq_poly_randtest(A, state, m, qbits);
        fmpq_poly_set_coeff_ui(A, 0, UWORD(0));
        arb_poly_randtest(b, state, 1 + n_randint(state, 20), rbits1, 5);

        fmpq_poly_exp_series(B, A, n);
        arb_poly_set_fmpq_poly(a, A, rbits1);
        arb_poly_exp_series(b, a, n, rbits2);

        if (!arb_poly_contains_fmpq_poly(b, B))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits2 = %wd\n", rbits2);

            flint_printf("A = "); fmpq_poly_print(A); flint_printf("\n\n");
            flint_printf("B = "); fmpq_poly_print(B); flint_printf("\n\n");

            flint_printf("a = "); arb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); arb_poly_printd(b, 15); flint_printf("\n\n");

            flint_abort();
        }

        fmpq_poly_clear(A);
        fmpq_poly_clear(B);

        arb_poly_clear(a);
        arb_poly_clear(b);
    }

    /* test aliasing */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n, qbits, rbits1, rbits2;
        fmpq_poly_t A;
        arb_poly_t a, b;

        qbits = 2 + n_randint(state, 100);
        rbits1 = 2 + n_randint(state, 100);
        rbits2 = 2 + n_randint(state, 100);

        if (n_randint(state, 100) == 0)
        {
            m = 1 + n_randint(state, 300);
            n = 1 + n_randint(state, 300);
        }
        else
        {
            m = 1 + n_randint(state, 20);
            n = 1 + n_randint(state, 20);
        }

        fmpq_poly_init(A);
        arb_poly_init(a);
        arb_poly_init(b);

        fmpq_poly_randtest(A, state, m, qbits);

        arb_poly_randtest(a, state, 1 + n_randint(state, 300), rbits1, 5);
        arb_poly_set_fmpq_poly(a, A, rbits1);

        arb_poly_exp_series(b, a, n, rbits2);
        arb_poly_exp_series(a, a, n, rbits2);

        if (!arb_poly_equal(a, b))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits2 = %wd\n", rbits2);

            flint_printf("A = "); fmpq_poly_print(A); flint_printf("\n\n");

            flint_printf("a = "); arb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); arb_poly_printd(b, 15); flint_printf("\n\n");

            flint_abort();
        }

        fmpq_poly_clear(A);
        arb_poly_clear(a);
        arb_poly_clear(b);
    }

    /* test that log(exp(f)) contains f */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n, qbits, rbits1, rbits2, rbits3;
        fmpq_poly_t A;
        arb_poly_t a, b, c;

        qbits = 2 + n_randint(state, 100);
        rbits1 = 2 + n_randint(state, 100);
        rbits2 = 2 + n_randint(state, 100);
        rbits3 = 2 + n_randint(state, 100);

        if (iter > 100)
        {
            arb_poly_newton_exp_cutoff = 5 + n_randint(state, 50);
        }

        if (n_randint(state, 100) == 0)
        {
            m = 1 + n_randint(state, 100);
            n = 1 + n_randint(state, 100);
        }
        else
        {
            m = 1 + n_randint(state, 20);
            n = 1 + n_randint(state, 20);
        }

        fmpq_poly_init(A);
        arb_poly_init(a);
        arb_poly_init(b);
        arb_poly_init(c);

        arb_poly_randtest(b, state, 1 + n_randint(state, 300), rbits1, 5);

        do {
            fmpq_poly_randtest(A, state, m, qbits);
            arb_poly_set_fmpq_poly(a, A, rbits1);
            arb_poly_exp_series(b, a, n, rbits3);
        } while (b->length == 0 || arb_contains_zero(b->coeffs));

        arb_poly_log_series(c, b, n, rbits2);

        fmpq_poly_truncate(A, n);

        if (!arb_poly_contains_fmpq_poly(c, A))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits2 = %wd\n", rbits2);
            flint_printf("bits3 = %wd\n", rbits3);

            flint_printf("A = "); fmpq_poly_print(A); flint_printf("\n\n");

            flint_printf("a = "); arb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); arb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); arb_poly_printd(c, 15); flint_printf("\n\n");

            flint_abort();
        }

        fmpq_poly_clear(A);
        arb_poly_clear(a);
        arb_poly_clear(b);
        arb_poly_clear(c);
    }

    TEST_FUNCTION_END(state);
}
