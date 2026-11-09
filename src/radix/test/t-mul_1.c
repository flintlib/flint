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
            flint_printf("radix %wu^%wd, n = %wd\n", DIGIT_RADIX(radix), radix->exp, n);
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

    /*
        Bases at the ends of the two-limb paths (B^2 does not fit a limb):
        B of FLINT_BITS / 2 + 1 bits, B below 2^(FLINT_BITS - 2) (the
        largest bases of the Barrett step), B of FLINT_BITS - 1 bits and
        B of FLINT_BITS bits, including powers of two; limbs and multipliers
        B - 1, B - 2, 1, 0 or random, longer operands
    */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        nn_ptr a, b, c;
        ulong B, x, cy;
        slong n, i;

        switch (n_randint(state, 10))
        {
            case 0: B = (UWORD(1) << (FLINT_BITS / 2)) + 1; break;
            case 1: B = (UWORD(1) << (FLINT_BITS / 2 + 1)) - 1; break;
            case 2: B = (UWORD(1) << (FLINT_BITS - 2)) - 1; break;
            case 3: B = (UWORD(1) << (FLINT_BITS - 3)) + 1; break;
            case 4: B = UWORD(1) << (FLINT_BITS - 2); break;
            case 5: B = (UWORD(1) << (FLINT_BITS - 1)) - 1; break;
            case 6: B = UWORD(1) << (FLINT_BITS - 1); break;
            case 7: B = (UWORD(1) << (FLINT_BITS - 1)) + 1; break;
            case 8: B = UWORD_MAX; break;
            default: B = n_randbits(state, FLINT_BITS / 2 + 1 + n_randint(state, FLINT_BITS / 2));
        }

        radix_init(radix, B, 1);

        n = 1 + n_randint(state, 40);

        a = flint_malloc(n * sizeof(ulong));
        b = flint_malloc(n * sizeof(ulong));
        c = flint_malloc((n + 1) * sizeof(ulong));

        for (i = 0; i < n; i++)
        {
            switch (n_randint(state, 5))
            {
                case 0: a[i] = B - 1; break;
                case 1: a[i] = B - 2; break;
                case 2: a[i] = n_randint(state, 2); break;
                default: a[i] = n_randint(state, B);
            }
        }

        switch (n_randint(state, 4))
        {
            case 0: x = B - 1; break;
            case 1: x = B - 2; break;
            case 2: x = 1; break;
            default: x = n_randint(state, B);
        }

        cy = radix_mul_1(b, a, n, x, radix);
        radix_mulmid_naive(c, a, n, &x, 1, 0, n + 1, radix);

        if (mpn_cmp(b, c, n) != 0 || cy != c[n])
        {
            flint_printf("FAIL: mul_1, two-limb bases\n");
            flint_printf("B = %wu, n = %wd, x = %wu\n", B, n, x);
            flint_printf("a = %{ulong*}\n", a, n);
            flint_printf("b = %{ulong*}\n", b, n);
            flint_printf("c = %{ulong*}\n", c, n + 1);
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
