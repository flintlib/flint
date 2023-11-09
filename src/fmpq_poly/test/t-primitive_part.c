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

TEST_FUNCTION_START(fmpq_poly_primitive_part, state)
{
    int i, result;
    ulong cflags = UWORD(0);

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g;
        fmpq_t x;

        fmpq_init(x);
        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_randtest(f, state, n_randint(state, 100), 200);
        fmpq_randtest(x, state, 100);

        fmpq_poly_scalar_mul_fmpq(f, f, x);

        fmpq_poly_primitive_part(g, f);
        fmpq_poly_primitive_part(f, f);

        cflags |= fmpq_poly_is_canonical(f) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(g) ? 0 : 2;
        result = (fmpq_poly_equal(f, g) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(f), flint_printf("\n\n");
            fmpq_poly_debug(g), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(x);
        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
    }

    /* Check that content(f) primitive_part(f) = +- f */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g;
        fmpq_t x, y;

        fmpq_init(x);
        fmpq_init(y);
        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_randtest(f, state, n_randint(state, 100), 200);
        fmpq_randtest(x, state, 100);

        fmpq_poly_scalar_mul_fmpq(f, f, x);

        fmpq_poly_content(y, f);
        fmpq_poly_primitive_part(g, f);
        fmpq_poly_scalar_mul_fmpq(g, g, y);

        if (!fmpq_poly_is_zero(f) && fmpz_sgn(f->coeffs + (f->length - 1)) < 0)
            fmpq_poly_neg(g, g);

        cflags |= fmpq_poly_is_canonical(f) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(g) ? 0 : 2;
        result = (fmpq_poly_equal(f, g) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(f), flint_printf("\n\n");
            fmpq_poly_debug(g), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
    }

    TEST_FUNCTION_END(state);
}
