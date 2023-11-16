/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "calcium.h"
#include "fexpr.h"

TEST_FUNCTION_START(fexpr_set_fmpz, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpz_t a, b;
        fexpr_t x;
        slong i;

        fmpz_init(a);
        fmpz_init(b);
        fexpr_init(x);

        for (i = 0; i < 5; i++)
        {
            fmpz_randtest(a, state, 300);
            fexpr_set_fmpz(x, a);
            fexpr_get_fmpz(b, x);

            if (!fmpz_equal(a, b))
            {
                flint_printf("FAIL\n\n");
                flint_printf("a = "); fmpz_print(a); printf("\n");
                flint_printf("b = "); fmpz_print(b); printf("\n");
                flint_printf("x = "); fexpr_print(x); printf("\n");
                flint_abort();
            }
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fexpr_clear(x);
    }

    TEST_FUNCTION_END(state);
}
