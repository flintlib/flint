/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_neg_uiui, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        ulong hi, lo;

        hi = n_randtest(state);
        lo = n_randtest(state);

        fmpz_init(a);
        fmpz_init(b);

        fmpz_set_ui(a, hi);
        fmpz_mul_2exp(a, a, FLINT_BITS);
        fmpz_add_ui(a, a, lo);
        fmpz_neg(a, a);

        fmpz_neg_uiui(b, hi, lo);

        result = fmpz_equal(a, b) && _fmpz_is_canonical(b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("hi = %wu\n", hi);
            flint_printf("lo = %wu\n", lo);
            flint_printf("a = "); fmpz_print(a); flint_printf("\n");
            flint_printf("b = "); fmpz_print(b); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
    }

    TEST_FUNCTION_END(state);
}
