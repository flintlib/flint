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

TEST_FUNCTION_START(fmpz_get_ui, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        ulong b, c;

        b = n_randtest(state);

        fmpz_init(a);

        fmpz_set_ui(a, b);
        c = fmpz_get_ui(a);

        result = (b == c) && _fmpz_is_canonical(a);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("b = %wd, c = %wd\n", b, c);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
    }

    TEST_FUNCTION_END(state);
}
