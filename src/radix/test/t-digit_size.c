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

TEST_FUNCTION_START(radix_digit_size, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        ulong c;
        ulong c1, c2;

        radix_init_randtest(radix, state);
        radix_randtest_limbs(&c, state, 1, radix);

        c1 = radix_size_digits_1(c, radix);

        c2 = 0;
        if (c != 0)
        {
            for (c2 = 1; c2 <= radix->exp; c2++)
            {
                if (n_pow(DIGIT_RADIX(radix), c2) > c)
                    break;
            }
        }

        if (c1 != c2)
        {
            flint_printf("FAIL: digit_size_1\n");
            flint_printf("radix %wu^%wd\n", DIGIT_RADIX(radix), radix->exp);
            flint_printf("c = %wu\n", c);
            flint_printf("c1 = %wu\n", c1);
            flint_printf("c2 = %wu\n", c2);
            flint_abort();
        }

        radix_clear(radix);
    }

    /* multi-limb numbers, possibly with zero top limbs */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        ulong x[6];
        ulong c1, c2;
        slong n, m, i;

        radix_init_randtest(radix, state);
        n = n_randint(state, 7);
        radix_randtest_limbs(x, state, n, radix);
        m = n_randint(state, n + 1);
        for (i = m; i < n; i++)
            x[i] = 0;

        c1 = radix_size_digits(x, n, radix);

        /* strip zero limbs and count the digits of the top limb */
        while (n > 0 && x[n - 1] == 0)
            n--;

        c2 = 0;
        if (n > 0)
        {
            for (c2 = 1; c2 <= radix->exp; c2++)
            {
                if (n_pow(DIGIT_RADIX(radix), c2) > x[n - 1])
                    break;
            }
            c2 += (n - 1) * radix->exp;
        }

        if (c1 != c2)
        {
            flint_printf("FAIL: size_digits\n");
            flint_printf("radix %wu^%wd\n", DIGIT_RADIX(radix), radix->exp);
            flint_printf("n = %wd\n", n);
            flint_printf("c1 = %wu\n", c1);
            flint_printf("c2 = %wu\n", c2);
            flint_abort();
        }

        radix_clear(radix);
    }

    TEST_FUNCTION_END(state);
}
