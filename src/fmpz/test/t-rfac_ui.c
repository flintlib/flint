/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_rfac_ui, state)
{
    int i, result;

    /* Check rf(x,a) * rf(x+a,b) = rf(x,a+b) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t x, xa, r1, r2, r1r2, r3;
        ulong a, b;

        fmpz_init(x);
        fmpz_init(xa);
        fmpz_init(r1);
        fmpz_init(r2);
        fmpz_init(r1r2);
        fmpz_init(r3);

        fmpz_randtest(x, state, 1 + n_randint(state, 200));
        a = n_randint(state, 100);
        b = n_randint(state, 100);
        fmpz_add_ui(xa, x, a);

        if (n_randint(state, 2))
        {
            fmpz_rfac_ui(r1, x, a);
        }
        else /* test aliasing */
        {
            fmpz_set(r1, x);
            fmpz_rfac_ui(r1, r1, a);
        }

        fmpz_rfac_ui(r2, xa, b);
        fmpz_rfac_ui(r3, x, a+b);

        fmpz_mul(r1r2, r1, r2);

        result = fmpz_equal(r1r2, r3) && _fmpz_is_canonical(r1) && _fmpz_is_canonical(r2) && _fmpz_is_canonical(r3);

        if (!result)
        {
            flint_printf("FAIL\n\n");
            flint_printf("x: "); fmpz_print(x); flint_printf("\n\n");
            flint_printf("a = %wu, b = %wu\n\n", a, b);
            flint_printf("rf(x,a): "); fmpz_print(r1); flint_printf("\n\n");
            flint_printf("rf(x+a,b): "); fmpz_print(r2); flint_printf("\n\n");
            flint_printf("rf(x,a+b): "); fmpz_print(r3); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(x);
        fmpz_clear(xa);
        fmpz_clear(r1);
        fmpz_clear(r2);
        fmpz_clear(r1r2);
        fmpz_clear(r3);
    }

    TEST_FUNCTION_END(state);
}
