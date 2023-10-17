/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

int sign(int a)
{
    return a > 0 ? 1 : a < 0 ? -1 : 0;
}

TEST_FUNCTION_START(fmpz_cmp2abs, state)
{
    slong i;

    for (i = 0; i < 200000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, b2;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(b2);

        fmpz_randtest(a, state, 400);
        fmpz_randtest(b, state, 400);
        fmpz_addmul_ui(a, b, 2);

        fmpz_mul_ui(b2, b, 2);
        if (sign(fmpz_cmp2abs(a, b)) != sign(fmpz_cmpabs(a, b2)))
        {
            flint_printf("FAIL i = %wd\n", i);
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(b2);
    }

    TEST_FUNCTION_END(state);
}
