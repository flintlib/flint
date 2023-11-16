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

TEST_FUNCTION_START(fmpzi_mul, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpzi_t x, y, z, s, t, u;

        fmpzi_init(x);
        fmpzi_init(y);
        fmpzi_init(z);
        fmpzi_init(s);
        fmpzi_init(t);
        fmpzi_init(u);

        fmpzi_randtest(x, state, 2000);
        fmpzi_randtest(y, state, 2000);
        fmpzi_randtest(z, state, 2000);
        fmpzi_randtest(s, state, 2000);
        fmpzi_randtest(t, state, 2000);

        /* check x * (y + z) = x * y + z * x */
        fmpzi_add(s, y, z);

        if (n_randint(state, 2))
            fmpzi_mul(s, x, s);
        else
            fmpzi_mul(s, s, x);

        fmpzi_mul(t, x, y);
        fmpzi_mul(u, z, x);
        fmpzi_add(t, t, u);

        if (!fmpzi_equal(s, t))
        {
            flint_printf("FAIL\n");
            flint_printf("x = "); fmpzi_print(x); printf("\n");
            flint_printf("y = "); fmpzi_print(y); printf("\n");
            flint_printf("z = "); fmpzi_print(z); printf("\n");
            flint_printf("s = "); fmpzi_print(s); printf("\n");
            flint_printf("t = "); fmpzi_print(t); printf("\n");
            flint_abort();
        }

        /* test squaring */
        fmpzi_randtest(x, state, 2000);
        fmpzi_randtest(y, state, 2000);
        fmpzi_randtest(s, state, 2000);
        fmpzi_randtest(t, state, 2000);

        fmpzi_set(y, x);
        fmpzi_mul(s, x, y);

        if (n_randint(state, 2))
        {
            fmpzi_mul(t, x, x);
        }
        else
        {
            fmpzi_set(t, x);
            fmpzi_mul(t, t, t);
        }

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
        fmpzi_clear(z);
        fmpzi_clear(s);
        fmpzi_clear(t);
        fmpzi_clear(u);
    }

    TEST_FUNCTION_END(state);
}
