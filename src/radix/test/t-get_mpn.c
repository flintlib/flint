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

TEST_FUNCTION_START(radix_get_mpn, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        nn_ptr a, b, c;
        slong n, l1, l2;

        radix_init_randtest(radix, state);
        n = n_randint(state, 200);

        a = flint_malloc(n * sizeof(ulong));
        b = flint_malloc(n * sizeof(ulong));
        c = flint_malloc(n * sizeof(ulong));

        radix_randtest_limbs(a, state, n, radix);

        l1 = radix_get_mpn_basecase(b, a, n, radix);
        l2 = radix_get_mpn(c, a, n, radix);

        if (l1 != l2 || mpn_cmp(b, c, l1) != 0)
        {
            flint_printf("FAIL\n");
            flint_printf("%{ulong*}\n", a, n);
            flint_printf("%{ulong*}\n", b, l1);
            flint_printf("%{ulong*}\n", c, l2);
            flint_abort();
        }

        radix_clear(radix);

        flint_free(a);
        flint_free(b);
        flint_free(c);
    }

    TEST_FUNCTION_END(state);
}
