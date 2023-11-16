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

TEST_FUNCTION_START(fmpzi_divexact, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpzi_t x, y, xy, q;
        int aliasing;
        slong bits;

        fmpzi_init(x);
        fmpzi_init(y);
        fmpzi_init(xy);
        fmpzi_init(q);

        if (n_randint(state, 10) == 0)
            bits = 2000;
        else
            bits = 200;

        fmpzi_randtest(x, state, 2 + n_randint(state, bits));
        do {
            fmpzi_randtest(y, state, 2 + n_randint(state, bits));
        } while (fmpzi_is_zero(y));

        fmpzi_randtest(q, state, 2 + n_randint(state, 200));

        fmpzi_mul(xy, x, y);

        aliasing = n_randint(state, 3);
        switch (aliasing)
        {
            case 0:
                fmpzi_divexact(q, xy, y);
                break;
            case 1:
                fmpzi_set(q, xy);
                fmpzi_divexact(q, q, y);
                break;
            case 2:
                fmpzi_set(q, y);
                fmpzi_divexact(q, xy, q);
                break;
        }

        if (!fmpzi_equal(q, x))
        {
            flint_printf("FAIL\n");
            flint_printf("aliasing = %d\n", aliasing);
            flint_printf("x = "); fmpzi_print(x); printf("\n");
            flint_printf("y = "); fmpzi_print(y); printf("\n");
            flint_printf("xy = "); fmpzi_print(xy); printf("\n");
            flint_printf("q = "); fmpzi_print(q); printf("\n");
            flint_abort();
        }

        fmpzi_clear(x);
        fmpzi_clear(y);
        fmpzi_clear(xy);
        fmpzi_clear(q);
    }

    TEST_FUNCTION_END(state);
}
