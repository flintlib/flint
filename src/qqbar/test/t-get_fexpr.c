/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fexpr.h"
#include "qqbar.h"

TEST_FUNCTION_START(qqbar_get_fexpr, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        qqbar_t x, y;
        fexpr_t e;

        qqbar_init(x);
        qqbar_init(y);
        fexpr_init(e);

        if (n_randint(state, 2))
            qqbar_randtest(x, state, 2, 10);
        else if (n_randint(state, 2))
            qqbar_randtest(x, state, 5, 100);
        else
            qqbar_randtest(x, state, 20, 20);

        qqbar_get_fexpr_repr(e, x);

        if (!qqbar_set_fexpr(y, e) || !qqbar_equal(x, y))
        {
            flint_printf("FAIL! (repr)\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("e = "); fexpr_print(e); flint_printf("\n\n");
            flint_abort();
        }

        qqbar_get_fexpr_root_indexed(e, x);
        qqbar_zero(y);

        if (!qqbar_set_fexpr(y, e) || !qqbar_equal(x, y))
        {
            flint_printf("FAIL! (root_indexed)\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("e = "); fexpr_print(e); flint_printf("\n\n");
            flint_abort();
        }

/*
        flint_printf("%wd\n", iter);
        flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
*/

        qqbar_get_fexpr_root_nearest(e, x);
        qqbar_zero(y);

/*
        flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
        flint_printf("e = "); fexpr_print(e); flint_printf("\n\n");
*/

        if (!qqbar_set_fexpr(y, e) || !qqbar_equal(x, y))
        {
            flint_printf("FAIL! (root_nearest)\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("e = "); fexpr_print(e); flint_printf("\n\n");
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(y);
        fexpr_clear(e);
    }

    TEST_FUNCTION_END(state);
}
