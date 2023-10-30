/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_neg_ui, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        ulong c;

        c = n_randtest(state);

        fmpz_init(a);
        fmpz_init(b);

        fmpz_set_ui(a, c);
        fmpz_neg(a, a);

        fmpz_neg_ui(b, c);

        result = fmpz_equal(a, b) && _fmpz_is_canonical(b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("c = %wu\n", c);
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
