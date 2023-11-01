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

TEST_FUNCTION_START(fmpz_poly_equal_fmpz, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b;
        fmpz_t c;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_init(c);

        fmpz_randtest(c, state, 1 + n_randint(state, 200));
        fmpz_poly_set_fmpz(b, c);

        fmpz_poly_randtest(a, state, n_randint(state, 6),
            1 + n_randint(state, 200));
        if (n_randint(state, 2))
            fmpz_poly_set_coeff_fmpz(a, 0, c);

        result = fmpz_poly_equal(a, b) == fmpz_poly_equal_fmpz(a, c);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fmpz_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_clear(c);
    }

    TEST_FUNCTION_END(state);
}
