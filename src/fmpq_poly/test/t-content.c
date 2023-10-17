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

TEST_FUNCTION_START(fmpq_poly_content, state)
{
    int i, result;

    /* Check that content(a f) = abs(a) content(f) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g;
        fmpq_t a, b, c;

        fmpq_poly_init(f);
        fmpq_poly_init(g);

        fmpq_init(a);
        fmpq_init(b);
        fmpq_init(c);

        fmpq_poly_randtest_not_zero(f, state, n_randint(state, 100) + 1, 100);
        fmpq_randtest_not_zero(a, state, 100);

        fmpq_poly_scalar_mul_fmpq(g, f, a);
        fmpq_poly_content(b, g);
        fmpq_poly_content(c, f);
        fmpq_mul(c, a, c);
        fmpq_abs(c, c);

        result = (fmpq_equal(b, c));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_print(f), flint_printf("\n\n");
            fmpq_poly_print(g), flint_printf("\n\n");
            fmpq_print(a), flint_printf("\n\n");
            fmpq_print(b), flint_printf("\n\n");
            fmpq_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_clear(a);
        fmpq_clear(b);
        fmpq_clear(c);
    }

    TEST_FUNCTION_END(state);
}
