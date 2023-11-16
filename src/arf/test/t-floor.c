/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arf.h"

TEST_FUNCTION_START(arf_floor, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, y;
        int result;

        arf_init(x);
        arf_init(y);

        arf_randtest_special(x, state, 2000, 100);
        arf_randtest_special(y, state, 2000, 100);

        arf_floor(y, x);

        result = 1;

        if (arf_is_int(x) || !arf_is_finite(x))
        {
            result = arf_equal(y, x);
        }
        else if (!arf_is_int(y))
        {
            result = 0;
        }
        else if (arf_cmp(y, x) >= 0)
        {
            result = 0;
        }
        else
        {
            arf_t s, t[3];

            /* check floor(x) - x + 1 > 0 */

            arf_init(s);
            arf_init(t[0]);
            arf_init(t[1]);
            arf_init(t[2]);

            arf_set(t[0], y);
            arf_neg(t[1], x);
            arf_one(t[2]);

            arf_sum(s, (arf_ptr) t, 3, 32, ARF_RND_DOWN);

            result = arf_sgn(s) > 0;

            arf_clear(s);
            arf_clear(t[0]);
            arf_clear(t[1]);
            arf_clear(t[2]);
        }

        if (!result)
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = "); arf_print(y); flint_printf("\n\n");
            flint_abort();
        }

        arf_floor(x, x);

        if (!arf_equal(x, y))
        {
            flint_printf("FAIL (aliasing)!\n");
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = "); arf_print(y); flint_printf("\n\n");
            flint_abort();
        }

        arf_clear(x);
        arf_clear(y);
    }

    TEST_FUNCTION_END(state);
}
