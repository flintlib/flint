/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_xgcd_partial, state)
{
    int i, result;

    /* Test co2*r1 - co1*r2 = r2_orig */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t co1, co2, f, g, t1, t2, L;

        fmpz_init(co1);
        fmpz_init(co2);
        fmpz_init(f);
        fmpz_init(g);
        fmpz_init(L);
        fmpz_init(t1);
        fmpz_init(t2);

        fmpz_randtest_unsigned(g, state, 200);
        fmpz_add_ui(g, g, 1);
        fmpz_randm(f, state, g);
        fmpz_randtest_unsigned(L, state, 200);

        fmpz_set(t2, g);
        fmpz_abs(t2, t2);

        fmpz_xgcd_partial(co2, co1, g, f, L);

        fmpz_mul(t1, co2, f);
        fmpz_submul(t1, co1, g);
        fmpz_abs(t1, t1);

        result = fmpz_equal(t1, t2);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("co1 = "), fmpz_print(co1), flint_printf("\n");
            flint_printf("co2 = "), fmpz_print(co2), flint_printf("\n");
            flint_printf("f = "), fmpz_print(f), flint_printf("\n");
            flint_printf("g = "), fmpz_print(g), flint_printf("\n");
            flint_printf("L = "), fmpz_print(L), flint_printf("\n");
            flint_printf("t1 = "), fmpz_print(t1), flint_printf("\n");
            flint_printf("t2 = "), fmpz_print(t2), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(co1);
        fmpz_clear(co2);
        fmpz_clear(f);
        fmpz_clear(g);
        fmpz_clear(L);
        fmpz_clear(t1);
        fmpz_clear(t2);
    }

    TEST_FUNCTION_END(state);
}
