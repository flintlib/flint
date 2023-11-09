/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"

TEST_FUNCTION_START(fmpq_div_fmpz, state)
{
    int i, result;

    /* Aliasing x = x/z */
    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, y;
        fmpz_t z;

        fmpq_init(x);
        fmpq_init(y);
        fmpz_init(z);

        fmpq_randtest(x, state, 200);
        while (fmpz_is_zero(z))
            fmpz_randtest(z, state, 200);

        fmpq_div_fmpz(y, x, z);
        fmpq_div_fmpz(x, x, z);

        result = (fmpq_is_canonical(x) && fmpq_is_canonical(y) && fmpq_equal(x, y));
        if (!result)
        {
            flint_printf("FAIL (alias):\n");
            flint_printf("x = "), fmpq_print(x), flint_printf("\n");
            flint_printf("y = "), fmpq_print(y), flint_printf("\n");
            flint_printf("z = "), fmpz_print(z), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(x);
        fmpq_clear(y);
        fmpz_clear(z);
    }

    /* Compare with fmpq_div */
    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, y, z;

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        fmpq_randtest(x, state, 200);
        while (fmpq_is_zero(z))
            fmpz_randtest(fmpq_numref(z), state, 200);

        fmpq_div_fmpz(y, x, fmpq_numref(z));
        fmpq_div(x, x, z);

        result = (fmpq_is_canonical(x) && fmpq_is_canonical(y) && fmpq_equal(x, y));
        if (!result)
        {
            flint_printf("FAIL (cmp):\n");
            flint_printf("x = "), fmpq_print(x), flint_printf("\n");
            flint_printf("y = "), fmpq_print(y), flint_printf("\n");
            flint_printf("z = "), fmpq_print(z), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    TEST_FUNCTION_END(state);
}
