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

TEST_FUNCTION_START(radix_n_divrem, state)
{
    slong iter, i;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        ulong d, x, q, r, q1, q2, q3, r1, r2, r3;
        n_div_precomp_t pre;

        d = n_randtest_not_zero(state);

        n_div_precomp_init(pre, d);

        for (i = 0; i < 10; i++)
        {
            x = n_randtest(state);

            q = x / d;
            r = x % d;

            q1 = n_divrem_precomp(&r1, x, d, pre);
            q2 = n_div_precomp(x, pre);
            r2 = n_rem_precomp(x, d, pre);

            if (pre->m == 0)
            {
                q3 = n_div_precomp_m0(x, pre);
                r3 = n_rem_precomp_m0(x, d, pre);
            }
            else if (pre->c == 0)
            {
                q3 = n_div_precomp_c0(x, pre);
                r3 = n_rem_precomp_c0(x, d, pre);
            }
            else if (x < UWORD_MAX && n_randint(state, 2))
            {
                q3 = n_div_precomp_c1_bounded(x, pre);
                r3 = n_rem_precomp_c1_bounded(x, d, pre);
            }
            else
            {
                q3 = n_div_precomp_c1(x, pre);
                r3 = n_rem_precomp_c1(x, d, pre);
            }

            if (q != q1 || q != q2 || q != q3 || r != r1 || r != r2 || r != r3)
            {
                flint_printf("FAIL: d = %wu, x = %wu\n", d, x);
                flint_printf("q = %wu, q1 = %wu, q2 = %wu, q3 = %wu\n", q, q1, q2, q3);
                flint_printf("r = %wu, r1 = %wu, r2 = %wu, r3 = %wu\n", r, r1, r2, r3);
                flint_abort();
            }
        }

    }

    TEST_FUNCTION_END(state);
}
