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

TEST_FUNCTION_START(radix_set_mpn, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        nn_ptr a, b, c, d;
        slong n, n2, n3, l1, l2, l3;

        radix_init_randtest(radix, state);

        do {
            n = n_randint(state, 100);
            n2 = radix_set_mpn_need_alloc(n, radix);
        } while (n2 > 500);
        n3 = n + 1;

        a = flint_malloc(n * sizeof(ulong));
        b = flint_malloc(n2 * sizeof(ulong));
        c = flint_malloc(n3 * sizeof(ulong));
        d = flint_malloc(n2 * sizeof(ulong));

        if (n != 0)
            flint_mpn_rrandom(a, state, n);

        radix_randtest_limbs(b, state, n2, radix);
        radix_randtest_limbs(d, state, n2, radix);

        if (n3 != 0)
            flint_mpn_rrandom(c, state, n3);

        l1 = radix_set_mpn_basecase(b, a, n, radix);
        l2 = radix_get_mpn(c, b, l1, radix);

        if (l2 > n || (l2 > 0 && mpn_cmp(a, c, l2) != 0) || (n > l2 && !mpn_zero_p(a + l2, n - l2)))
        {
            flint_printf("FAIL set_mpn\n");
            flint_printf("%{ulong*}\n", a, n);
            flint_printf("%{ulong*}\n", b, l1);
            flint_printf("%{ulong*}\n", c, l2);
            flint_abort();
        }

        l3 = radix_set_mpn(d, a, n, radix);

        if (l3 != l1 || (l1 > 0 && mpn_cmp(b, d, l1) != 0))
        {
            flint_printf("FAIL set_mpn\n");
            flint_printf("%{ulong*}\n", a, n);
            flint_printf("set_mpn_basecase:\n");
            flint_printf("%{ulong*}\n", b, l1);
            flint_printf("set_mpn:\n");
            flint_printf("%{ulong*}\n", d, l3);
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
