/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpzi.h"

TEST_FUNCTION_START(fmpzi_add_sub, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        /* check x - y = x + (-y) */
        fmpzi_t x, y, s, t;

        fmpzi_init(x);
        fmpzi_init(y);
        fmpzi_init(s);
        fmpzi_init(t);

        fmpzi_randtest(x, state, 100);
        fmpzi_randtest(y, state, 100);
        fmpzi_randtest(s, state, 100);
        fmpzi_randtest(t, state, 100);

        fmpzi_sub(s, x, y);
        fmpzi_neg(t, y);
        fmpzi_add(t, x, t);

        if (!fmpzi_equal(s, t))
        {
            flint_printf("FAIL\n");
            flint_printf("x = "); fmpzi_print(x); printf("\n");
            flint_printf("y = "); fmpzi_print(y); printf("\n");
            flint_printf("s = "); fmpzi_print(s); printf("\n");
            flint_printf("t = "); fmpzi_print(t); printf("\n");
            flint_abort();
        }

        fmpzi_clear(x);
        fmpzi_clear(y);
        fmpzi_clear(s);
        fmpzi_clear(t);
    }

    TEST_FUNCTION_END(state);
}
