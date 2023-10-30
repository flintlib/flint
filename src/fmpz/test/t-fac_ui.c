/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_fac_ui, state)
{
    slong i, n;
    fmpz_t x;
    fmpz_t y;

    fmpz_init(x);
    fmpz_init(y);

    /* Twice to check demotion */
    for (n = 0; n < 2; n++)
    {
        fmpz_set_ui(y, UWORD(1));

        for (i = 0; i < 100; i++)
        {
            fmpz_fac_ui(x, i);
            fmpz_mul_ui(y, y, FLINT_MAX(1, i));
            if (!fmpz_equal(x, y) || !_fmpz_is_canonical(x))
            {
                flint_printf("FAIL: %wd\n", i);
                fmpz_print(x);
                flint_printf("\n");
                fmpz_print(y);
                flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }
    }

    fmpz_clear(x);
    fmpz_clear(y);

    TEST_FUNCTION_END(state);
}
