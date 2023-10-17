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

TEST_FUNCTION_START(fmpz_poly_zero_coeffs, state)
{
    int i, result;

    /* Check that zeroing [0,len/2) and [len/2,len) sets a to zero */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a;
        slong len;

        fmpz_poly_init(a);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 200);
        len = a->length;

        fmpz_poly_zero_coeffs(a, -23, len/2);
        fmpz_poly_zero_coeffs(a, len/2, len + 42);

        result = (fmpz_poly_length(a) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
    }

    TEST_FUNCTION_END(state);
}
