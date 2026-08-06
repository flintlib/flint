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

/* whether a[0] is a square unit modulo the digit radix p */
static int
_is_square_unit2(nn_srcptr a, const radix_t radix)
{
    ulong p = DIGIT_RADIX(radix);
    ulong r = a[0] % p;

    if (r == 0)
        return 0;
    if (p == 2)
        return (a[0] & 7) == 1;
    return n_sqrtmod(r, p) != 0;
}

TEST_FUNCTION_START(radix_sqrtmod_bn, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        nn_ptr a, s, c;
        slong an, n, i;

        radix_init_randtest_prime(radix, state);
        if (DIGIT_RADIX(radix) == 2 && radix->exp <= 2)
        {
            radix_clear(radix);
            radix_init(radix, 2, 3);
        }

        an = 1 + n_randint(state, 50);
        n = 1 + n_randint(state, 50);

        a = flint_malloc(an * sizeof(ulong));
        s = flint_malloc(n * sizeof(ulong));
        c = flint_malloc(n * sizeof(ulong));

        radix_randtest_limbs(a, state, an, radix);

        if (radix_sqrtmod_bn(s, a, an, n, radix))
        {
            /* check s^2 == a (mod B^n) */
            radix_mulmid(c, s, n, s, n, 0, n, radix);

            for (i = 0; i < n; i++)
            {
                ulong ai = (i < an) ? a[i] : 0;
                if (c[i] != ai)
                {
                    flint_printf("FAIL: sqrtmod_bn (s^2 != a)\n");
                    flint_printf("radix %wu ^ %u = %wu\n", DIGIT_RADIX(radix), radix->exp, LIMB_RADIX(radix));
                    flint_printf("%{ulong*}\n", a, an);
                    flint_printf("%{ulong*}\n", s, n);
                    flint_printf("%{ulong*}\n", c, n);
                    flint_abort();
                }
            }
        }
        else if (_is_square_unit2(a, radix))
        {
            flint_printf("FAIL: sqrtmod_bn (should be a square)\n");
            flint_printf("radix %wu ^ %u = %wu\n", DIGIT_RADIX(radix), radix->exp, LIMB_RADIX(radix));
            flint_printf("%{ulong*}\n", a, an);
            flint_abort();
        }

        radix_clear(radix);

        flint_free(a);
        flint_free(s);
        flint_free(c);
    }

    TEST_FUNCTION_END(state);
}
