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

TEST_FUNCTION_START(radix_mul_1, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        nn_ptr a, b, c;
        ulong x;
        slong n;
        ulong cy;

        radix_init_randtest(radix, state);

        n = 1 + n_randint(state, 5);

        a = flint_malloc(n * sizeof(ulong));
        b = flint_malloc(n * sizeof(ulong));
        c = flint_malloc((n + 1) * sizeof(ulong));

        radix_randtest_limbs(a, state, n, radix);
        radix_randtest_limbs(b, state, n, radix);
        radix_randtest_limbs(c, state, n + 1, radix);
        radix_randtest_limbs(&x, state, 1, radix);

        if (n_randint(state, 2))
        {
            cy = radix_mul_1(b, a, n, x, radix);
        }
        else
        {
            flint_mpn_copyi(b, a, n);
            cy = radix_mul_1(b, b, n, x, radix);
        }

        radix_mulmid_naive(c, a, n, &x, 1, 0, n + 1, radix);

        if (mpn_cmp(b, c, n) != 0 || cy != c[n])
        {
            flint_printf("FAIL: mul_1\n");
            flint_printf("radix %wu^%wd, n = %wd\n", radix->B.n, radix->exp, n);
            flint_printf("a = %{ulong*}\n", a, n);
            flint_printf("b = %{ulong*}\n", b, n);
            flint_printf("c = %{ulong*}\n", c, n + 1);
            flint_printf("x = %wu\n", x);
            flint_printf("cy = %wu\n", cy);
            flint_abort();
        }

        radix_clear(radix);

        flint_free(a);
        flint_free(b);
        flint_free(c);
    }

    TEST_FUNCTION_END(state);
}
