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

TEST_FUNCTION_START(fmpz_poly_primitive_part, state)
{
    int i, result;

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_randtest(g, state, n_randint(state, 100), 200);

        fmpz_poly_primitive_part(f, g);
        fmpz_poly_primitive_part(g, g);

        result = (fmpz_poly_equal(f, g));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(f), flint_printf("\n\n");
            fmpz_poly_print(g), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
    }

    /* Check that content(f) primitive_part(f) = sgn(lead(f)) f */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g;
        fmpz_t c;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_init(c);
        fmpz_poly_randtest_not_zero(f, state, n_randint(state, 100) + 1, 200);

        fmpz_poly_content(c, f);
        if (fmpz_sgn(f->coeffs + f->length - 1) < 0)
            fmpz_neg(c, c);
        fmpz_poly_primitive_part(g, f);
        fmpz_poly_scalar_mul_fmpz(g, g, c);

        result = (fmpz_poly_equal(f, g));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(f), flint_printf("\n\n");
            fmpz_poly_print(g), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(c);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
    }

    TEST_FUNCTION_END(state);
}
