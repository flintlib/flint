/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2009 William Hart

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

TEST_FUNCTION_START(fmpz_poly_pow_binomial, state)
{
    int i, result;

    /* Check aliasing of a and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b;
        ulong exp;

        fmpz_poly_init(a);
        fmpz_poly_init2(b, 2);

        fmpz_randtest(b->coeffs, state, 100);
        fmpz_randtest_not_zero(b->coeffs + 1, state, 100);
        _fmpz_poly_set_length(b, 2);

        exp = n_randtest(state) % UWORD(100);

        fmpz_poly_pow_binomial(a, b, exp);
        fmpz_poly_pow_binomial(b, b, exp);

        result = (fmpz_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL(1):\n");
            flint_printf("exp = %wu\n", exp);
            flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
    }

    /* Compare with fmpz_poly_pow */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b;
        ulong exp;

        fmpz_poly_init(a);
        fmpz_poly_init2(b, 2);

        fmpz_randtest(b->coeffs, state, 100);
        fmpz_randtest_not_zero(b->coeffs + 1, state, 100);
        _fmpz_poly_set_length(b, 2);

        exp = n_randtest(state) % UWORD(100);

        fmpz_poly_pow_binomial(a, b, exp);
        fmpz_poly_pow(b, b, exp);

        result = (fmpz_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL(2):\n");
            flint_printf("exp = %wu\n", exp);
            flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
    }

    TEST_FUNCTION_END(state);
}
