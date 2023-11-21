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

TEST_FUNCTION_START(acb_poly_compose_series, state)
{
    slong iter;

    for (iter = 0; iter < 3000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong qbits1, qbits2, rbits1, rbits2, rbits3, n;
        fmpq_poly_t A, B, C;
        acb_poly_t a, b, c, d;

        qbits1 = 2 + n_randint(state, 200);
        qbits2 = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        rbits3 = 2 + n_randint(state, 200);
        n = 2 + n_randint(state, 25);

        fmpq_poly_init(A);
        fmpq_poly_init(B);
        fmpq_poly_init(C);

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);
        acb_poly_init(d);

        fmpq_poly_randtest(A, state, 1 + n_randint(state, 25), qbits1);
        fmpq_poly_randtest(B, state, 1 + n_randint(state, 25), qbits2);
        fmpq_poly_set_coeff_ui(B, 0, 0);
        fmpq_poly_compose_series(C, A, B, n);

        acb_poly_randtest(c, state, 1 + n_randint(state, 20), rbits1, 4);

        acb_poly_set_fmpq_poly(a, A, rbits1);
        acb_poly_set_fmpq_poly(b, B, rbits2);
        acb_poly_compose_series(c, a, b, n, rbits3);

        if (!acb_poly_contains_fmpq_poly(c, C))
        {
            flint_printf("FAIL\n\n");
            flint_printf("n = %wd, bits3 = %wd\n", n, rbits3);

            flint_printf("A = "); fmpq_poly_print(A); flint_printf("\n\n");
            flint_printf("B = "); fmpq_poly_print(B); flint_printf("\n\n");
            flint_printf("C = "); fmpq_poly_print(C); flint_printf("\n\n");

            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_poly_printd(c, 15); flint_printf("\n\n");

            flint_abort();
        }

        acb_poly_set(d, a);
        acb_poly_compose_series(d, d, b, n, rbits3);
        if (!acb_poly_equal(d, c))
        {
            flint_printf("FAIL (aliasing 1)\n\n");
            flint_abort();
        }

        acb_poly_set(d, b);
        acb_poly_compose_series(d, a, d, n, rbits3);
        if (!acb_poly_equal(d, c))
        {
            flint_printf("FAIL (aliasing 2)\n\n");
            flint_abort();
        }

        if (a->length > 0 && b->length > 0)
        {
            int any_finite;
            slong i, k, l, m;

            /* randomize coefficients to set to indeterminate value */
            /* if the value is equal to the length we don't set an indeterminate value */
            k = n_randint(state, a->length + 1);
            l = 1 + n_randint(state, b->length);

            if (k < a->length)
              acb_indeterminate(a->coeffs + k);
            else
              k = n; // k doesn't affect number of finite coefficients

            if (l < b->length)
              acb_indeterminate(b->coeffs + l);
            else
              l = n; // l doesn't affect number of finite coefficients

            acb_poly_compose_series(d, a, b, n, rbits3);

            /* up to this all coefficients should be finite */
            m = FLINT_MIN(FLINT_MIN(k, l), d->length);

            /* check that coefficients after m are all non-finite */
            any_finite = 0;
            for (i = m; i < d->length; i++)
                any_finite |= acb_is_finite(d->coeffs + i);

            if (any_finite)
            {
                flint_printf("FAIL (non-finite 1)\n\n");
                flint_abort();
            }

            /* check that coefficients up to m are all finite and
               contain the expected result */
            if (!_acb_vec_is_finite(d->coeffs, m))
            {
                flint_printf("FAIL (non-finite 2)\n\n");
                flint_abort();
            }

            fmpq_poly_truncate(C, m);
            acb_poly_truncate(d, m);

            if (!acb_poly_contains_fmpq_poly(d, C))
            {
                flint_printf("FAIL (non-finite 3)\n\n");
                flint_abort();
            }
        }

        fmpq_poly_clear(A);
        fmpq_poly_clear(B);
        fmpq_poly_clear(C);

        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(c);
        acb_poly_clear(d);
    }

    TEST_FUNCTION_END(state);
}
