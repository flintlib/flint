/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "radix.h"

TEST_FUNCTION_START(radix_div_exp, state)
{
    slong iter, i;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        ulong x, q, r, e;

        radix_init_randtest(radix, state);
        e = radix->exp;

        for (i = 0; i < 100; i++)
        {
            switch (n_randint(state, 4))
            {
                case 0: x = n_randint(state, 1000); break;
                case 1: x = n_randtest(state); break;
                case 2: x = (UWORD(1) << (FLINT_BITS - 2)) - 1 - n_randint(state, 100); break;
                default: x = UWORD_MAX - n_randint(state, 100); break;
            }

            q = radix_div_exp(x, radix);
            if (q != x / e || radix_divrem_exp(&r, x, radix) != q || r != x % e)
            {
                flint_printf("FAIL: radix_div_exp\n");
                flint_printf("exp = %wu, x = %wu, q = %wu\n", e, x, q);
                flint_abort();
            }
        }

        radix_clear(radix);
    }

    TEST_FUNCTION_END(state);
}
