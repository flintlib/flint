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

TEST_FUNCTION_START(qqbar_guess, state)
{
    slong iter;

    for (iter = 0; iter < 200 * 0.1 * flint_test_multiplier(); iter++)
    {
        qqbar_t x, y;
        acb_t z;
        slong prec, max_deg, max_bits;

        qqbar_init(x);
        qqbar_init(y);
        acb_init(z);

        if (n_randint(state, 2))
            qqbar_randtest(x, state, 4, 100);
        else
            qqbar_randtest(x, state, 12, 10);

        max_deg = qqbar_degree(x) + n_randint(state, 2);
        max_bits = qqbar_height_bits(x) + n_randint(state, 100);

        for (prec = 64; ; prec *= 2)
        {
            qqbar_get_acb(z, x, prec);

            if (qqbar_guess(y, z, max_deg, max_bits, 0, prec) && qqbar_equal(x, y))
                break;

            if (prec > 10000)
            {
                flint_printf("FAIL!\n");
                flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
                flint_printf("prec = %wd\n\n", prec);
                flint_abort();
            }
        }

        qqbar_clear(x);
        qqbar_clear(y);
        acb_clear(z);
    }

    TEST_FUNCTION_END(state);
}
