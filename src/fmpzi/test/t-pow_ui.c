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

TEST_FUNCTION_START(fmpzi_pow_ui, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpzi_t x, y, s, t, u, v;
        ulong e, f;

        fmpzi_init(x);
        fmpzi_init(y);
        fmpzi_init(s);
        fmpzi_init(t);
        fmpzi_init(u);
        fmpzi_init(v);

        fmpzi_randtest(x, state, 100);
        fmpzi_randtest(y, state, 100);
        fmpzi_randtest(s, state, 100);
        fmpzi_randtest(t, state, 100);
        fmpzi_randtest(u, state, 100);

        e = n_randint(state, 30);
        f = n_randint(state, 30);

        if (n_randint(state, 2))
        {
            fmpzi_pow_ui(s, x, e);
        }
        else
        {
            fmpzi_set(s, x);
            fmpzi_pow_ui(s, s, e);
        }

        fmpzi_pow_ui(t, x, f);
        fmpzi_mul(v, s, t);

        fmpzi_pow_ui(u, x, e + f);

        if (!fmpzi_equal(v, u))
        {
            flint_printf("FAIL\n");
            flint_printf("x = "); fmpzi_print(x); printf("\n");
            flint_printf("y = "); fmpzi_print(y); printf("\n");
            flint_printf("s = "); fmpzi_print(s); printf("\n");
            flint_printf("t = "); fmpzi_print(t); printf("\n");
            flint_printf("u = "); fmpzi_print(u); printf("\n");
            flint_printf("v = "); fmpzi_print(v); printf("\n");
            flint_abort();
        }

        fmpzi_clear(x);
        fmpzi_clear(y);
        fmpzi_clear(s);
        fmpzi_clear(t);
        fmpzi_clear(u);
        fmpzi_clear(v);
    }

    TEST_FUNCTION_END(state);
}
