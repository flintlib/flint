/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"

TEST_FUNCTION_START(fmpq_height, state)
{
    int i;

    for (i = 0; i < 10000; i++)
    {
        fmpq_t x;
        fmpz_t h;
        flint_bitcnt_t b;

        fmpz_init(h);
        fmpq_init(x);
        fmpq_randtest(x, state, 200);

        fmpq_height(h, x);
        b = fmpq_height_bits(x);

        if (fmpz_bits(h) != b)
        {
            flint_printf("FAIL!\n");
            flint_printf("x: ");
            fmpq_print(x);
            flint_printf("\nh: ");
            fmpz_print(h);
            flint_printf("\nb: %wd\n", b);
        }

        fmpq_clear(x);
        fmpz_clear(h);
    }

    TEST_FUNCTION_END(state);
}
