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
_is_square_unit(nn_srcptr a, const radix_t radix)
{
    ulong p = DIGIT_RADIX(radix);
    ulong r = a[0] % p;

    if (r == 0)
        return 0;                       /* not a unit */
    if (p == 2)
        return (a[0] & 7) == 1;         /* odd squares are 1 mod 8 */
    return n_sqrtmod(r, p) != 0;        /* quadratic residue mod p */
}

TEST_FUNCTION_START(radix_rsqrtmod_bn, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        nn_ptr a, b, c, d;
        slong an, n;

        radix_init_randtest_prime(radix, state);
        if (DIGIT_RADIX(radix) == 2 && radix->exp <= 2)
        {
            radix_clear(radix);
            radix_init(radix, 2, 3);
        }

        an = 1 + n_randint(state, 50);
        n = 1 + n_randint(state, 50);

        a = flint_malloc(an * sizeof(ulong));
        b = flint_malloc(n * sizeof(ulong));
        c = flint_malloc(n * sizeof(ulong));
        d = flint_malloc(n * sizeof(ulong));

        radix_randtest_limbs(a, state, an, radix);

        if (radix_rsqrtmod_bn(b, a, an, n, radix))
        {
            /* check a * b^2 == 1 (mod B^n) */
            radix_mulmid(c, b, n, b, n, 0, n, radix);
            radix_mulmid(d, a, FLINT_MIN(an, n), c, n, 0, n, radix);

            if (d[0] != 1 || !flint_mpn_zero_p(d + 1, n - 1))
            {
                flint_printf("FAIL: rsqrtmod_bn (a b^2 != 1)\n");
                flint_printf("radix %wu ^ %u = %wu\n", DIGIT_RADIX(radix), radix->exp, LIMB_RADIX(radix));
                flint_printf("%{ulong*}\n", a, an);
                flint_printf("%{ulong*}\n", b, n);
                flint_printf("%{ulong*}\n", d, n);
                flint_abort();
            }
        }
        else if (_is_square_unit(a, radix))
        {
            flint_printf("FAIL: rsqrtmod_bn (should be a square)\n");
            flint_printf("radix %wu ^ %u = %wu\n", DIGIT_RADIX(radix), radix->exp, LIMB_RADIX(radix));
            flint_printf("%{ulong*}\n", a, an);
            flint_abort();
        }

        radix_clear(radix);

        flint_free(a);
        flint_free(b);
        flint_free(c);
        flint_free(d);
    }

    TEST_FUNCTION_END(state);
}
