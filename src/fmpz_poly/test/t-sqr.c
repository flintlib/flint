/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_poly_sqr, state)
{
    int i, result;

    /* Check aliasing of a and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(a, state, n_randint(state, 50), 500);

        fmpz_poly_set(b, a);
        fmpz_poly_sqr(c, b);
        fmpz_poly_sqr(b, b);

        result = (fmpz_poly_equal(b, c));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fmpz_poly_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Compare with mul */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 500);

        fmpz_poly_sqr(b, a);
        fmpz_poly_mul(c, a, a);

        result = (fmpz_poly_equal(b, c));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fmpz_poly_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Check _fmpz_poly_sqr directly */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        slong len;
        fmpz_poly_t a, out1, out2;

        len = n_randint(state, 100) + 1;
        fmpz_poly_init(a);
        fmpz_poly_init(out1);
        fmpz_poly_init(out2);
        fmpz_poly_randtest(a, state, len, 200);

        fmpz_poly_sqr(out1, a);
        fmpz_poly_fit_length(a, a->alloc + n_randint(state, 10));
        a->length = a->alloc;
        fmpz_poly_fit_length(out2, 2 * a->length - 1);
        _fmpz_poly_sqr(out2->coeffs, a->coeffs, a->length);
        _fmpz_poly_set_length(out2, 2 * a->length - 1);
        _fmpz_poly_normalise(out2);

        result = (fmpz_poly_equal(out1, out2));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(out1), flint_printf("\n\n");
            fmpz_poly_print(out2), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(out1);
        fmpz_poly_clear(out2);
    }

    TEST_FUNCTION_END(state);
}
