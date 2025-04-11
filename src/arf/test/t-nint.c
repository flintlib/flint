/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arf.h"

TEST_FUNCTION_START(arf_nint, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        arf_t x, y, s, t;

        arf_init(x);
        arf_init(y);
        arf_init(s);
        arf_init(t);

        arf_randtest_special(x, state, 200, 100);
        arf_randtest_special(y, state, 200, 100);

        if (n_randint(state, 2))
        {
            arf_nint(y, x);
        }
        else
        {
            arf_set(y, x);
            arf_nint(y, y);
        }

        if (arf_cmpabs_2exp_si(x, -3) < 0)
            arf_zero(s);
        else if (arf_is_int(x) || arf_is_special(x))
            arf_floor(s, x);
        else
        {
            /* nint(x) = floor(x+0.5) - isint((2*x-1)/4) */
            arf_set_d(s, 0.5);
            arf_add(s, s, x, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_floor(s, s);
            arf_mul_2exp_si(t, x, 1);
            arf_sub_ui(t, t, 1, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_mul_2exp_si(t, t, -2);
            arf_sub_ui(s, s, arf_is_int(t), ARF_PREC_EXACT, ARF_RND_DOWN);
        }

        if (!arf_equal(y, s))
        {
            flint_printf("FAIL\n");
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = "); arf_print(y); flint_printf("\n\n");
            flint_printf("s = "); arf_print(s); flint_printf("\n\n");
            flint_abort();
        }

        arf_clear(x);
        arf_clear(y);
        arf_clear(s);
        arf_clear(t);
    }

    TEST_FUNCTION_END(state);
}

