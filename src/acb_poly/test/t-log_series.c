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

TEST_FUNCTION_START(acb_poly_log_series, state)
{
    slong iter;

    /* compare with fmpq_poly */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n, qbits, rbits1, rbits2;
        fmpq_poly_t A, B;
        acb_poly_t a, b;

        qbits = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 20);
        n = 1 + n_randint(state, 20);

        fmpq_poly_init(A);
        fmpq_poly_init(B);

        acb_poly_init(a);
        acb_poly_init(b);

        fmpq_poly_randtest_not_zero(A, state, m, qbits);
        fmpq_poly_set_coeff_ui(A, 0, UWORD(1));

        fmpq_poly_log_series(B, A, n);
        acb_poly_set_fmpq_poly(a, A, rbits1);
        acb_poly_randtest(b, state, 1 + n_randint(state, 20), rbits1, 5);
        acb_poly_log_series(b, a, n, rbits2);

        if (!acb_poly_contains_fmpq_poly(b, B))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits2 = %wd\n", rbits2);

            flint_printf("A = "); fmpq_poly_print(A); flint_printf("\n\n");
            flint_printf("B = "); fmpq_poly_print(B); flint_printf("\n\n");

            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");

            flint_abort();
        }

        fmpq_poly_clear(A);
        fmpq_poly_clear(B);

        acb_poly_clear(a);
        acb_poly_clear(b);
    }

    /* test aliasing */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n, qbits, rbits1, rbits2;
        fmpq_poly_t A;
        acb_poly_t a, b;

        qbits = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 20);
        n = 1 + n_randint(state, 20);

        fmpq_poly_init(A);
        acb_poly_init(a);
        acb_poly_init(b);

        do {
            fmpq_poly_randtest_not_zero(A, state, m, qbits);
        } while (fmpz_sgn(A->coeffs + 0) <= 0);

        acb_poly_set_fmpq_poly(a, A, rbits1);

        acb_poly_log_series(b, a, n, rbits2);
        acb_poly_log_series(a, a, n, rbits2);

        if (!acb_poly_equal(a, b))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits2 = %wd\n", rbits2);

            flint_printf("A = "); fmpq_poly_print(A); flint_printf("\n\n");

            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");

            flint_abort();
        }

        fmpq_poly_clear(A);
        acb_poly_clear(a);
        acb_poly_clear(b);
    }

    /* test that exp(log(f)) contains f */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n, qbits, rbits1, rbits2, rbits3;
        fmpq_poly_t A;
        acb_poly_t a, b, c;

        qbits = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        rbits3 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 20);
        n = 1 + n_randint(state, 20);

        fmpq_poly_init(A);
        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);

        do {
            fmpq_poly_randtest_not_zero(A, state, m, qbits);
        } while (fmpz_sgn(A->coeffs + 0) <= 0);

        acb_poly_set_fmpq_poly(a, A, rbits1);

        acb_poly_randtest(b, state, 1 + n_randint(state, 20), rbits1, 5);
        acb_poly_log_series(b, a, n, rbits2);
        acb_poly_exp_series_basecase(c, b, n, rbits3);

        fmpq_poly_truncate(A, n);

        if (!acb_poly_contains_fmpq_poly(c, A))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits2 = %wd\n", rbits2);
            flint_printf("bits3 = %wd\n", rbits3);

            flint_printf("A = "); fmpq_poly_print(A); flint_printf("\n\n");

            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_poly_printd(c, 15); flint_printf("\n\n");

            flint_abort();
        }

        fmpq_poly_clear(A);
        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(c);
    }

    TEST_FUNCTION_END(state);
}
