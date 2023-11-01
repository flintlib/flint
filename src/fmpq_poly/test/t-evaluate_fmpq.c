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
#include "fmpq.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_evaluate_fmpq, state)
{
    int i, result;

    /* Check aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        fmpq_t x, y;
        fmpq_poly_t f;

        fmpz_init(a);
        fmpz_init(b);
        fmpq_init(x);
        fmpq_init(y);
        fmpq_poly_init(f);
        fmpq_poly_randtest(f, state, n_randint(state, 80), 100);
        fmpz_randtest(a, state, 80);
        fmpz_randtest_not_zero(b, state, 80);
        fmpz_set(fmpq_numref(x), a);
        fmpz_set(fmpq_denref(x), a);
        fmpq_canonicalise(x);

        fmpq_poly_evaluate_fmpq(y, f, x);
        fmpq_poly_evaluate_fmpq(x, f, x);

        result = (fmpq_equal(x, y));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_print(a), flint_printf("\n\n");
            fmpz_print(b), flint_printf("\n\n");
            fmpq_poly_debug(f), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_poly_clear(f);
    }

    /* Check that (f+g)(a) = f(a) + g(a) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        fmpq_t x, y, z;
        fmpq_poly_t f, g;

        fmpz_init(a);
        fmpz_init(b);
        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);
        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_randtest(f, state, n_randint(state, 80), 100);
        fmpq_poly_randtest(g, state, n_randint(state, 80), 100);
        fmpz_randtest(a, state, 80);
        fmpz_randtest_not_zero(b, state, 80);
        fmpz_set(fmpq_numref(x), a);
        fmpz_set(fmpq_denref(x), a);
        fmpq_canonicalise(x);

        fmpq_poly_evaluate_fmpq(y, f, x);
        fmpq_poly_evaluate_fmpq(z, g, x);
        fmpq_add(y, y, z);
        fmpq_poly_add(f, f, g);
        fmpq_poly_evaluate_fmpq(z, f, x);

        result = (fmpq_equal(y, z));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_print(a), flint_printf("\n\n");
            fmpz_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
    }

    /* Check that (f*g)(a) = f(a) * g(a) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        fmpq_t x, y, z;
        fmpq_poly_t f, g;

        fmpz_init(a);
        fmpz_init(b);
        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);
        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_randtest(f, state, n_randint(state, 50), 80);
        fmpq_poly_randtest(g, state, n_randint(state, 50), 80);
        fmpz_randtest(a, state, 80);
        fmpz_randtest_not_zero(b, state, 80);
        fmpz_set(fmpq_numref(x), a);
        fmpz_set(fmpq_denref(x), a);
        fmpq_canonicalise(x);

        fmpq_poly_evaluate_fmpq(y, f, x);
        fmpq_poly_evaluate_fmpq(z, g, x);
        fmpq_mul(y, y, z);
        fmpq_poly_mul(f, f, g);
        fmpq_poly_evaluate_fmpq(z, f, x);

        result = (fmpq_equal(y, z));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_print(a), flint_printf("\n\n");
            fmpz_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
    }

    TEST_FUNCTION_END(state);
}
