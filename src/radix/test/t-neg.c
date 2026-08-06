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

TEST_FUNCTION_START(radix_neg, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        nn_ptr a, b, c, d;
        slong an;
        ulong cy1, cy2;

        radix_init_randtest(radix, state);

        an = 1 + n_randint(state, 5);

        a = flint_malloc(an * sizeof(ulong));
        b = flint_malloc(an * sizeof(ulong));
        c = flint_malloc(an * sizeof(ulong));
        d = flint_malloc(an * sizeof(ulong));

        radix_randtest_limbs(a, state, an, radix);
        radix_randtest_limbs(c, state, an, radix);

        flint_mpn_zero(b, an);

        cy1 = radix_neg(c, a, an, radix);
        cy2 = radix_sub(d, b, an, a, an, radix);

        if (cy1 != cy2 || mpn_cmp(c, d, an) != 0)
        {
            flint_printf("FAIL: neg\n");
            flint_printf("%{ulong*}\n", a, an);
            flint_printf("%{ulong*}\n", b, an);
            flint_printf("%{ulong*}\n", c, an);
            flint_printf("%{ulong*}\n", d, an);
            flint_printf("%wu\n%wu\n\n", cy1, cy2);
            flint_abort();
        }

        flint_mpn_copyi(c, a, an);
        cy1 = radix_neg(c, c, an, radix);

        if (cy1 != cy2 || mpn_cmp(c, d, an) != 0)
        {
            flint_printf("FAIL: neg (aliasing)\n");
            flint_printf("%{ulong*}\n", a, an);
            flint_printf("%{ulong*}\n", b, an);
            flint_printf("%{ulong*}\n", c, an);
            flint_printf("%{ulong*}\n", d, an);
            flint_printf("%wu\n%wu\n\n", cy1, cy2);
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
