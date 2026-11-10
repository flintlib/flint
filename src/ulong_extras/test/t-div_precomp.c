/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_div_precomp, state)
{
    slong iter, i;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        ulong d, x, q, r, q1, q2, q3, r1, r2, r3, m;
        unsigned int e, c;
        n_div_precomp_t pre, pre2;

        d = n_randtest_not_zero(state);

        n_div_precomp_init(pre, d);

        /* the data against its definition, by a two-word division:
           m = floor(2^(FLINT_BITS + e) / d) + 1 - c */
        e = FLINT_BIT_COUNT(d) - 1;
        if ((d & (d - 1)) == 0)
        {
            m = 0;
            c = 0;
        }
        else
        {
            udiv_qrnnd(q, r, UWORD(1) << e, UWORD(0), d);
            c = (d - r >= (UWORD(1) << e));
            m = q + 1 - c;
        }
        if (pre->m != m || pre->e != e || pre->c != c)
        {
            flint_printf("FAIL: n_div_precomp_init, d = %wu\n", d);
            flint_printf("m = %wu, %wu; e = %u, %u; c = %u, %u\n",
                         pre->m, m, pre->e, e, pre->c, c);
            flint_abort();
        }

        /* the same from the inverse of d (n_div_precomp_init may not use it) */
        n_div_precomp_init_preinv(pre2, d, n_preinvert_limb(d));
        if (pre2->m != m || pre2->e != e || pre2->c != c)
        {
            flint_printf("FAIL: n_div_precomp_init_preinv, d = %wu\n", d);
            flint_printf("m = %wu, %wu; e = %u, %u; c = %u, %u\n",
                         pre2->m, m, pre2->e, e, pre2->c, c);
            flint_abort();
        }

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
