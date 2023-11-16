/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "qqbar.h"

TEST_FUNCTION_START(qqbar_root_of_unity, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        qqbar_t x, y;
        slong p, p2;
        ulong q, q2;
        int ok;

        qqbar_init(x);
        qqbar_init(y);

        q = 1 + n_randint(state, 30);
        p = n_randint(state, 1000);
        p -= 500;

        qqbar_root_of_unity(x, p, q);
        ok = qqbar_is_root_of_unity(&p2, &q2, x);
        qqbar_root_of_unity(y, p2, q2);

        if (!ok || !qqbar_equal(x, y) || n_gcd(FLINT_ABS(p2), q2) != 1 || !(0 <= p2 && p2 < (slong) 2 * q2))
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("p, p2 = %wd %wd\n\n", p, p2);
            flint_printf("q, q2 = %wu %wu\n\n", q, q2);
            flint_abort();
        }

        if (q <= 10)
        {
            qqbar_pow_ui(y, x, q);

            if (!qqbar_is_one(y))
            {
                flint_printf("FAIL!\n");
                flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
                flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
                flint_printf("p, p2 = %wd %wd\n\n", p, p2);
                flint_printf("q, q2 = %wu %wu\n\n", q, q2);
                flint_abort();
            }
        }

        qqbar_clear(x);
        qqbar_clear(y);
    }

    TEST_FUNCTION_END(state);
}
