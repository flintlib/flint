/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "double_extras.h"
#include "qqbar.h"

TEST_FUNCTION_START(qqbar_set_re_im_d, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        double x, y;
        qqbar_t z;
        acb_t a, b;
        int ok;

        qqbar_init(z);
        acb_init(a);
        acb_init(b);

        x = d_randtest_special(state, -1100, 1100);
        y = d_randtest_special(state, -1100, 1100);

        ok =  qqbar_set_re_im_d(z, x, y);

        qqbar_get_acb(a, z, 53);

        acb_set_d_d(b, x, y);

        if (ok)
        {
            if (!acb_equal(a, b))
            {
                flint_printf("FAIL!\n");
                flint_printf("z = "); qqbar_print(z); flint_printf("\n\n");
                flint_printf("a = "); acb_print(a); flint_printf("\n\n");
                flint_printf("b = "); acb_print(b); flint_printf("\n\n");
                flint_abort();
            }
        }
        else
        {
            if (acb_is_finite(b))
            {
                flint_printf("FAIL!\n");
                flint_printf("z = "); qqbar_print(z); flint_printf("\n\n");
                flint_printf("a = "); acb_print(a); flint_printf("\n\n");
                flint_printf("b = "); acb_print(b); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(a);
        acb_clear(b);
        qqbar_clear(z);
    }

    TEST_FUNCTION_END(state);
}
