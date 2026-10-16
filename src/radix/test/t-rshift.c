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

TEST_FUNCTION_START(radix_rshift, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        nn_ptr a, b, c;
        slong n;
        ulong cy1, cy2;
        unsigned int e;
        ulong be;

        radix_init_randtest(radix, state);
        while (radix->exp == 1)
        {
            radix_clear(radix);
            radix_init_randtest(radix, state);
        }

        n = 1 + n_randint(state, 5);
        e = 1 + n_randint(state, radix->exp - 1);

        a = flint_malloc(n * sizeof(ulong));
        b = flint_malloc(n * sizeof(ulong));
        c = flint_malloc(n * sizeof(ulong));

        radix_randtest_limbs(a, state, n, radix);
        radix_randtest_limbs(b, state, n, radix);
        radix_randtest_limbs(c, state, n, radix);

        if (n_randint(state, 2))
        {
            cy1 = radix_rshift_digits(b, a, n, e, radix);
        }
        else
        {
            flint_mpn_copyi(b, a, n);
            cy1 = radix_rshift_digits(b, b, n, e, radix);
        }

        be = n_pow(DIGIT_RADIX(radix), e);
        cy2 = radix_divrem_1(c, a, n, be, radix) * n_pow(DIGIT_RADIX(radix), radix->exp - e);

        if (mpn_cmp(b, c, n) != 0 || cy1 != cy2)
        {
            flint_printf("FAIL: rshift\n");
            flint_printf("radix %wu^%wd, n = %wd\n", DIGIT_RADIX(radix), radix->exp, n);
            flint_printf("a = %{ulong*}\n", a, n);
            flint_printf("b = %{ulong*}\n", b, n);
            flint_printf("c = %{ulong*}\n", c, n);
            flint_printf("be = %wu\n", be);
            flint_printf("cy1 = %wu, cy2 = %wu\n", cy2);
            flint_abort();
        }

        radix_clear(radix);

        flint_free(a);
        flint_free(b);
        flint_free(c);
    }

    TEST_FUNCTION_END(state);
}
