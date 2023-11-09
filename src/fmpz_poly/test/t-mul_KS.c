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
#include "fmpz_poly.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_poly_mul_KS, state)
{
    int i, result;

    /* Check aliasing of a and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, state, n_randint(state, 50), 200);
        fmpz_poly_randtest(c, state, n_randint(state, 50), 200);
        fmpz_poly_mul_KS(a, b, c);
        fmpz_poly_mul_KS(b, b, c);

        result = (fmpz_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, state, n_randint(state, 50), 200);
        fmpz_poly_randtest(c, state, n_randint(state, 50), 200);

        fmpz_poly_mul_KS(a, b, c);
        fmpz_poly_mul_KS(c, b, c);

        result = (fmpz_poly_equal(a, c));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Check aliasing of b and c */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, state, n_randint(state, 50), 200);
        fmpz_poly_set(c, b);

        fmpz_poly_mul_KS(a, b, b);
        fmpz_poly_mul_KS(c, b, c);

        result = (fmpz_poly_equal(a, c));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Compare with mul_classical */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c, d;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(d);
        fmpz_poly_randtest(b, state, n_randint(state, 50), 200);
        fmpz_poly_randtest(c, state, n_randint(state, 50), 200);

        fmpz_poly_mul_KS(a, b, c);
        fmpz_poly_mul_classical(d, b, c);

        result = (fmpz_poly_equal(a, d));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(d), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(d);
    }

    /* Compare with mul_classical unsigned */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c, d;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(d);
        fmpz_poly_randtest_unsigned(b, state, n_randint(state, 50), 200);
        fmpz_poly_randtest_unsigned(c, state, n_randint(state, 50), 200);

        fmpz_poly_mul_KS(a, b, c);
        fmpz_poly_mul_classical(d, b, c);

        result = (fmpz_poly_equal(a, d));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(d), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(d);
    }

    /* Check _fmpz_poly_mul_KS directly */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        slong len1, len2;
        fmpz_poly_t a, b, out1, out2;

        len1 = n_randint(state, 100) + 1;
        len2 = n_randint(state, 100) + 1;
        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(out1);
        fmpz_poly_init(out2);
        fmpz_poly_randtest(a, state, len1, 200);
        fmpz_poly_randtest(b, state, len2, 200);

        fmpz_poly_mul_KS(out1, a, b);
        fmpz_poly_fit_length(a, a->alloc + n_randint(state, 10));
        fmpz_poly_fit_length(b, b->alloc + n_randint(state, 10));
        a->length = a->alloc;
        b->length = b->alloc;
        fmpz_poly_fit_length(out2, a->length + b->length - 1);
        _fmpz_poly_mul_KS(out2->coeffs, a->coeffs, a->length,
                                        b->coeffs, b->length);
        _fmpz_poly_set_length(out2, a->length + b->length - 1);
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
        fmpz_poly_clear(b);
        fmpz_poly_clear(out1);
        fmpz_poly_clear(out2);
    }

    TEST_FUNCTION_END(state);
}
