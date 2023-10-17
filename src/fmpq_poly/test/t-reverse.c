/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_reverse, state)
{
    int i, result;

    /* Aliasing */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;
        slong n;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(b, state, n_randint(state, 100), 200);
        n = n_randint(state, 150);

        fmpq_poly_reverse(a, b, n);
        fmpq_poly_reverse(b, b, n);

        result = (fmpq_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("n = %wd\n", n);
            flint_printf("a = "), fmpq_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpq_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    /* Correctness (?) */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;
        slong j, len, n;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(b, state, n_randint(state, 100), 200);
        n = n_randint(state, 150);

        len = FLINT_MIN(n, b->length);
        if (len)
        {
            fmpq_poly_fit_length(a, n);
            for (j = 0; j < len; j++)
                fmpz_set(a->coeffs + (n - len) + j, b->coeffs + (len - 1 - j));
            fmpz_set(a->den, b->den);
            _fmpq_poly_set_length(a, n);
            fmpq_poly_canonicalise(a);
        }

        fmpq_poly_reverse(b, b, n);

        result = (fmpq_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("n = %wd\n", n);
            flint_printf("a = "), fmpq_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpq_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    TEST_FUNCTION_END(state);
}
