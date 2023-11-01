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

TEST_FUNCTION_START(fmpq_one, state)
{
    int i, result;

    /* x == 1 * x */
    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, y, z;

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        fmpq_randtest(x, state, 200);
        fmpq_one(y);

        fmpq_mul(z, y, x);

        result = fmpq_is_canonical(z) && fmpq_equal(x, z);
        if (!result)
        {
            flint_printf("FAIL:\n");
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

    /* x/x == 1 */
    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, y;

        fmpq_init(x);
        fmpq_init(y);

        while (fmpq_is_zero(x))
            fmpq_randtest(x, state, 200);

        fmpq_div(y, x, x);

        result = fmpq_is_canonical(y) && fmpq_is_one(y);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("x = "), fmpq_print(x), flint_printf("\n");
            flint_printf("y = "), fmpq_print(y), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(x);
        fmpq_clear(y);
    }

    TEST_FUNCTION_END(state);
}
