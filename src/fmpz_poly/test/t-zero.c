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

TEST_FUNCTION_START(fmpz_poly_zero, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a;

        fmpz_poly_init(a);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 100);
        fmpz_poly_zero(a);

        result = (fmpz_poly_is_zero(a));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
    }

    TEST_FUNCTION_END(state);
}
