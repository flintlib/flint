/*
    Copyright (C) 1509 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_poly_evaluate_divconquer_fmpq, state)
{
    int i, result;

    /* Check aliasing */
    for (i = 0; i < 300 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        fmpq_t x, y;
        fmpz_poly_t f;

        fmpz_init(a);
        fmpz_init(b);
        fmpq_init(x);
        fmpq_init(y);
        fmpz_poly_init(f);
        fmpz_poly_randtest(f, state, n_randint(state, 100), 150);
        fmpz_randtest(a, state, 100);
        fmpz_randtest_not_zero(b, state, 100);
        fmpz_set(fmpq_numref(x), a);
        fmpz_set(fmpq_denref(x), b);
        fmpq_canonicalise(x);

        fmpz_poly_evaluate_divconquer_fmpq(y, f, x);
        fmpz_poly_evaluate_divconquer_fmpq(x, f, x);

        result = (fmpq_equal(x, y));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_print(a), flint_printf("\n\n");
            fmpz_print(b), flint_printf("\n\n");
            fmpz_poly_print(f), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpq_clear(x);
        fmpq_clear(y);
        fmpz_poly_clear(f);
    }

    /* Check that (f+g)(a) = f(a) + g(a) */
    for (i = 0; i < 300 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        fmpq_t x, y, z;
        fmpz_poly_t f, g;

        fmpz_init(a);
        fmpz_init(b);
        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_randtest(f, state, n_randint(state, 100), 150);
        fmpz_poly_randtest(g, state, n_randint(state, 100), 150);
        fmpz_randtest(a, state, 100);
        fmpz_randtest_not_zero(b, state, 100);
        fmpz_set(fmpq_numref(x), a);
        fmpz_set(fmpq_denref(x), b);
        fmpq_canonicalise(x);

        fmpz_poly_evaluate_divconquer_fmpq(y, f, x);
        fmpz_poly_evaluate_divconquer_fmpq(z, g, x);
        fmpq_add(y, y, z);
        fmpz_poly_add(f, f, g);
        fmpz_poly_evaluate_divconquer_fmpq(z, f, x);

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
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
    }

    /* Check that (f*g)(a) = f(a) * g(a) */
    for (i = 0; i < 300 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        fmpq_t x, y, z;
        fmpz_poly_t f, g;

        fmpz_init(a);
        fmpz_init(b);
        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_randtest(f, state, n_randint(state, 50), 100);
        fmpz_poly_randtest(g, state, n_randint(state, 50), 100);
        fmpz_randtest(a, state, 100);
        fmpz_randtest_not_zero(b, state, 100);
        fmpz_set(fmpq_numref(x), a);
        fmpz_set(fmpq_denref(x), b);
        fmpq_canonicalise(x);

        fmpz_poly_evaluate_divconquer_fmpq(y, f, x);
        fmpz_poly_evaluate_divconquer_fmpq(z, g, x);
        fmpq_mul(y, y, z);
        fmpz_poly_mul(f, f, g);
        fmpz_poly_evaluate_divconquer_fmpq(z, f, x);

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
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
    }

    TEST_FUNCTION_END(state);
}
