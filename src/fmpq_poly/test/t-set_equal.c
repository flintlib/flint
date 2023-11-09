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
#include "fmpq.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_set_equal, state)
{
    int i, result;

    /* equal polynomials */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);

        fmpq_poly_set(b, a);

        result = (fmpq_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpq_poly_debug(a), flint_printf("\n");
            flint_printf("alloc = %wd\nlength = %wd\n\n", a->alloc, a->length);
            flint_printf("b = "), fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("alloc = %wd\nlength = %wd\n\n", b->alloc, b->length);
            flint_printf("equal(a, b) = %d\n", result);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;
        slong coeff = (slong) n_randint(state, 100);
        fmpq_t x1, x2;
        fmpz_t x1fmpz;

        fmpq_init(x1);
        fmpq_init(x2);
        fmpz_init(x1fmpz);
        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);

        fmpq_poly_set(b, a);

        fmpq_poly_get_coeff_fmpq(x2, b, coeff);
        do
        {
            fmpz_randtest(x1fmpz, state, 200);
            fmpz_set(fmpq_numref(x1), x1fmpz);
            fmpz_set_si(fmpq_denref(x1), 1);
        } while (fmpq_equal(x1, x2));
        fmpq_poly_set_coeff_fmpq(b, coeff, x1);

        result = (!fmpq_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpq_poly_debug(a), flint_printf("\n");
            flint_printf("alloc = %wd\nlength = %wd\n\n", a->alloc, a->length);
            flint_printf("b = "), fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("alloc = %wd\nlength = %wd\n\n", b->alloc, b->length);
            flint_printf("!equal(a, b) = %d\n", result);
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(x1);
        fmpq_clear(x2);
        fmpz_clear(x1fmpz);
        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    TEST_FUNCTION_END(state);
}
