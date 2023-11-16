/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "qqbar.h"

TEST_FUNCTION_START(qqbar_sgn_re, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        qqbar_t x, y, z;
        int s1, s2, s3, s4;

        qqbar_init(x);
        qqbar_init(y);
        qqbar_init(z);

        qqbar_randtest(x, state, 4, 10);
        qqbar_i(y);
        qqbar_mul(y, y, x);

        s1 = qqbar_sgn_re(x);
        s2 = qqbar_sgn_im(y);

        if (s1 != s2)
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("s1 = %d", s1); flint_printf("\n\n");
            flint_printf("s2 = %d", s2); flint_printf("\n\n");
            flint_abort();
        }

        qqbar_randtest(x, state, 4, 10);
        qqbar_randtest(y, state, 4, 10);
        if (n_randint(state, 2))
            qqbar_add(x, x, y);
        else
            qqbar_mul(x, x, y);

        qqbar_i(y);
        qqbar_mul(y, y, x);

        s1 = qqbar_sgn_re(x);
        s2 = qqbar_sgn_im(y);

        if (s1 != s2)
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("s1 = %d", s1); flint_printf("\n\n");
            flint_printf("s2 = %d", s2); flint_printf("\n\n");
            flint_abort();
        }

        qqbar_randtest_real(x, state, 2, 100);
        qqbar_randtest_real(y, state, 2, 100);

        s1 = qqbar_sgn_re(x);
        s2 = qqbar_sgn_re(y);

        qqbar_i(z);
        qqbar_mul(y, y, z);
        qqbar_add(x, x, y);

        s3 = qqbar_sgn_re(x);
        s4 = qqbar_sgn_im(x);

        if (s1 != s3 || s2 != s4)
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("s1 = %d", s1); flint_printf("\n\n");
            flint_printf("s2 = %d", s2); flint_printf("\n\n");
            flint_printf("s3 = %d", s3); flint_printf("\n\n");
            flint_printf("s4 = %d", s4); flint_printf("\n\n");
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(y);
        qqbar_clear(z);
    }

    TEST_FUNCTION_END(state);
}
