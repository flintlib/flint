/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"
#include "arb_poly.h"

/* hack: avoid overflow since exp currently uses mpfr */
void
fmpq_poly_randtest_small(fmpq_poly_t A, flint_rand_t state, slong len, slong bits)
{
    fmpq_poly_randtest(A, state, len, bits);
    if (A->length > 0)
    {
        bits = _fmpz_vec_max_bits(A->coeffs, A->length);
        bits = FLINT_ABS(bits);
        fmpz_mul_2exp(A->den, A->den, bits);
        _fmpq_poly_normalise(A);
    }
}

TEST_FUNCTION_START(arb_poly_exp_series_basecase, state)
{
    slong iter;

    /* compare with fmpq_poly */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n, qbits, rbits1, rbits2;
        fmpq_poly_t A, B;
        arb_poly_t a, b;

        qbits = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 40);
        n = 1 + n_randint(state, 40);

        fmpq_poly_init(A);
        fmpq_poly_init(B);

        arb_poly_init(a);
        arb_poly_init(b);

        fmpq_poly_randtest(A, state, m, qbits);
        fmpq_poly_set_coeff_ui(A, 0, UWORD(0));

        fmpq_poly_exp_series(B, A, n);
        arb_poly_set_fmpq_poly(a, A, rbits1);
        arb_poly_exp_series_basecase(b, a, n, rbits2);

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

        qbits = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 40);
        n = 1 + n_randint(state, 40);

        fmpq_poly_init(A);
        arb_poly_init(a);
        arb_poly_init(b);

        fmpq_poly_randtest_small(A, state, m, qbits);

        arb_poly_set_fmpq_poly(a, A, rbits1);

        arb_poly_exp_series_basecase(b, a, n, rbits2);
        arb_poly_exp_series_basecase(a, a, n, rbits2);

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

        qbits = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        rbits3 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 40);
        n = 1 + n_randint(state, 40);

        fmpq_poly_init(A);
        arb_poly_init(a);
        arb_poly_init(b);
        arb_poly_init(c);

        do {
            fmpq_poly_randtest_small(A, state, m, qbits);
            arb_poly_set_fmpq_poly(a, A, rbits1);
            arb_poly_exp_series_basecase(b, a, n, rbits3);
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
