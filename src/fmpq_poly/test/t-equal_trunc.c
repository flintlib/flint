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

TEST_FUNCTION_START(fmpq_poly_equal_trunc, state)
{
    int i, result;

    /* equal polynomials */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;
        fmpq_t c;
        slong n, j;

        fmpq_init(c);
        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpq_poly_randtest(b, state, n_randint(state, 100), 200);
        n = n_randint(state, 100);

        for (j = 0; j < n; j++)
        {
           fmpq_poly_get_coeff_fmpq(c, a, j);
           fmpq_poly_set_coeff_fmpq(b, j, c);
        }

        result = (fmpq_poly_equal_trunc(a, b, n));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpq_poly_debug(a), flint_printf("\n");
            flint_printf("alloc = %wd\nlength = %wd\n\n", a->alloc, a->length);
            flint_printf("b = "), fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("alloc = %wd\nlength = %wd\n\n", b->alloc, b->length);
            flint_printf("equal(a, b) = %d\n", result);
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(c);
        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    /* unequal polynomials */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;
        fmpq_t c;
        slong n, m, j;

        fmpq_init(c);
        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpq_poly_randtest(b, state, n_randint(state, 100), 200);
        n = n_randint(state, 100) + 1;

        for (j = 0; j < n; j++)
        {
           fmpq_poly_get_coeff_fmpq(c, a, j);
           fmpq_poly_set_coeff_fmpq(b, j, c);
        }

        m = n_randint(state, n);
        fmpq_poly_get_coeff_fmpq(c, a, m);
        fmpq_add_si(c, c, 1);
        fmpq_poly_set_coeff_fmpq(b, m, c);

        result = (!fmpq_poly_equal_trunc(a, b, n));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpq_poly_debug(a), flint_printf("\n");
            flint_printf("alloc = %wd\nlength = %wd\n\n", a->alloc, a->length);
            flint_printf("b = "), fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("alloc = %wd\nlength = %wd\n\n", b->alloc, b->length);
            flint_printf("equal(a, b) = %d\n", result);
            flint_printf("n = %wd, m = %wd\n", n, m);
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(c);
        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    TEST_FUNCTION_END(state);
}
