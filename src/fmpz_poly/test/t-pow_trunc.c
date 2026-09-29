/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_poly_pow_trunc, state)
{
    int i, result;

    /* Check aliasing of a and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b;
        slong n;
        ulong exp;

        n   = n_randtest(state) % 10;
        exp = n_randtest(state) % 100;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_randtest(b, state, n_randint(state, 100), n);

        fmpz_poly_pow_trunc(a, b, exp, n);
        fmpz_poly_pow_trunc(b, b, exp, n);

        result = (fmpz_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("n = %wd\n", n);
            flint_printf("exp = %wu\n", exp);
            flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
    }

    /* Compare with powering followed truncating */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;
        slong n;
        ulong exp;

        n   = n_randtest(state) % 10;
        exp = n_randtest(state) % 50;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, state, n_randint(state, 50), n);

        fmpz_poly_pow(a, b, exp);
        fmpz_poly_truncate(a, n);
        fmpz_poly_pow_trunc(c, b, exp, n);

        result = (fmpz_poly_equal(a, c));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("n   = %wd\n", n);
            flint_printf("exp = %wu\n", exp);
            flint_printf("b = "), fmpz_poly_print(b), flint_printf("\n\n");
            flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n\n");
            flint_printf("c = "), fmpz_poly_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Check the underscore method against powering followed by
       truncation, including e = 1 and zero-padded input */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;
        slong n, len;
        ulong exp;

        n   = 1 + n_randint(state, 20);
        exp = 1 + n_randint(state, 10);
        len = n_randint(state, n + 1);

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, state, len, 1 + n_randint(state, 100));

        fmpz_poly_pow(a, b, exp);
        fmpz_poly_truncate(a, n);

        fmpz_poly_fit_length(b, n);
        _fmpz_vec_zero(b->coeffs + b->length, n - b->length);
        fmpz_poly_fit_length(c, n);
        _fmpz_poly_pow_trunc(c->coeffs, b->coeffs, exp, n);
        _fmpz_poly_set_length(c, n);
        _fmpz_poly_normalise(c);

        result = (fmpz_poly_equal(a, c));
        if (!result)
        {
            flint_printf("FAIL (underscore):\n");
            flint_printf("n   = %wd\n", n);
            flint_printf("exp = %wu\n", exp);
            flint_printf("b = "), fmpz_poly_print(b), flint_printf("\n\n");
            flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n\n");
            flint_printf("c = "), fmpz_poly_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    TEST_FUNCTION_END(state);
}
