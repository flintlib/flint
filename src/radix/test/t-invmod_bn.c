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

TEST_FUNCTION_START(radix_invmod_bn, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        nn_ptr a, b, c;
        slong an, n;

        radix_init_randtest(radix, state);

        an = 1 + n_randint(state, 10);
        n = 1 + n_randint(state, 10);

        a = flint_malloc(an * sizeof(ulong));
        b = flint_malloc(n * sizeof(ulong));
        c = flint_malloc(n * sizeof(ulong));

        radix_randtest_limbs(a, state, an, radix);
        radix_randtest_limbs(b, state, n, radix);

        if (radix_invmod_bn(b, a, an, n, radix))
        {
            radix_mulmid(c, b, n, a, FLINT_MIN(an, n), 0, n, radix);

            if (c[0] != 1 || !flint_mpn_zero_p(c + 1, n - 1))
            {
                flint_printf("FAIL: invmod_bn\n");
                flint_printf("radix %wu ^ %u = %wu\n", DIGIT_RADIX(radix), radix->exp, LIMB_RADIX(radix));
                flint_printf("%{ulong*}\n", a, an);
                flint_printf("%{ulong*}\n", b, n);
                flint_printf("%{ulong*}\n", c, n);
                flint_abort();
            }
        }
        else
        {
            if (n_gcd(a[0], LIMB_RADIX(radix)) == 1)
            {
                flint_printf("FAIL: invmod_bn (should be invertible)\n");
                flint_printf("radix %wu ^ %u = %wu\n", DIGIT_RADIX(radix), radix->exp, LIMB_RADIX(radix));
                flint_printf("%{ulong*}\n", a, an);
            }
        }

        radix_clear(radix);

        flint_free(a);
        flint_free(b);
        flint_free(c);
    }

    TEST_FUNCTION_END(state);
}
