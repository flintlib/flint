/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_asinh_series, state)
{
    int i, result;
    ulong cflags = UWORD(0);

    /* Check aliasing of a and c */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;
        slong n = n_randint(state, 50) + 1;

        fmpq_poly_init(a);
        fmpq_poly_init(b);

        fmpq_poly_randtest_not_zero(a, state, n_randint(state, 50) + 1, 50);
        fmpq_poly_set_coeff_ui(a, 0, UWORD(0));

        fmpq_poly_canonicalise(a);

        fmpq_poly_asinh_series(b, a, n);
        fmpq_poly_asinh_series(a, a, n);

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(b) ? 0 : 2;
        result = (fmpq_poly_equal(a, b) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    /* Check asinh(A) = atanh(A/sqrt(1+A^2)) */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t A, B, asinhA, atanhB;
        slong n = n_randint(state, 80) + 1;

        fmpq_poly_init(A);
        fmpq_poly_init(B);
        fmpq_poly_init(asinhA);
        fmpq_poly_init(atanhB);

        fmpq_poly_randtest_not_zero(A, state, n_randint(state, 80) + 1, 80);
        fmpq_poly_set_coeff_ui(A, 0, UWORD(0));

        fmpq_poly_mullow(B, A, A, n);
        fmpq_poly_set_coeff_ui(B, 0, UWORD(1));
        fmpq_poly_invsqrt_series(B, B, n);
        fmpq_poly_mullow(B, A, B, n);

        fmpq_poly_asinh_series(asinhA, A, n);
        fmpq_poly_atanh_series(atanhB, B, n);

        cflags |= fmpq_poly_is_canonical(asinhA) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(atanhB) ? 0 : 2;
        result = (fmpq_poly_equal(asinhA, atanhB) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("A = "), fmpq_poly_debug(A), flint_printf("\n\n");
            flint_printf("B = "), fmpq_poly_debug(B), flint_printf("\n\n");
            flint_printf("asinh(A) = "), fmpq_poly_debug(asinhA), flint_printf("\n\n");
            flint_printf("atanh(B) = "), fmpq_poly_debug(atanhB), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(A);
        fmpq_poly_clear(B);
        fmpq_poly_clear(asinhA);
        fmpq_poly_clear(atanhB);
    }

    TEST_FUNCTION_END(state);
}
