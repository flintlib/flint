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

TEST_FUNCTION_START(radix_sqrt, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        nn_ptr a, s, s0;
        slong an, sn, sn0;
        int exact, want_exact;

        radix_t radix;
        radix_init_randtest(radix, state);

        sn0 = 1 + n_randint(state, n_randint(state, 100) ? 20 : 150);

        a = flint_malloc(sizeof(ulong) * (2 * sn0 + 1));
        s0 = flint_malloc(sizeof(ulong) * sn0);

        radix_randtest_limbs(s0, state, sn0, radix);
        if (s0[sn0 - 1] == 0)
            s0[sn0 - 1] = 1 + n_randint(state, LIMB_RADIX(radix) - 1);

        /* a = s0^2 + d with 0 <= d <= 2*s0, so floor(sqrt(a)) = s0 */
        radix_mul(a, s0, sn0, s0, sn0, radix);
        an = 2 * sn0;
        MPN_NORM(a, an);

        want_exact = n_randint(state, 2);
        if (!want_exact)
        {
            /* add d in [1, 2*s0]: either a small d or 2*s0 itself */
            a[an] = 0;
            if (n_randint(state, 2))
            {
                ulong dmax = LIMB_RADIX(radix) - 1;
                if (sn0 == 1 && 2 * s0[0] < dmax)
                    dmax = 2 * s0[0];
                ulong d = 1 + n_randint(state, dmax);
                radix_add(a, a, an + 1, &d, 1, radix);
            }
            else
            {
                nn_ptr t = flint_malloc(sizeof(ulong) * (sn0 + 1));
                t[sn0] = radix_add(t, s0, sn0, s0, sn0, radix);
                radix_add(a, a, an + 1, t, sn0 + 1, radix);
                flint_free(t);
            }
            an = an + 1;
            MPN_NORM(a, an);
        }

        sn = (an + 1) / 2;
        s = flint_malloc(sizeof(ulong) * sn);

        exact = radix_sqrt(s, a, an, radix);

        if (exact != want_exact)
        {
            flint_printf("FAIL: wrong exactness flag: radix %wu^%wd, an = %wd, exact = %d\n", DIGIT_RADIX(radix), radix->exp, an, exact);
            flint_printf("a = %{ulong*}\n", a, an);
            flint_abort();
        }

        if (sn != sn0 || mpn_cmp(s, s0, sn0) != 0)
        {
            flint_printf("FAIL: wrong root: radix %wu^%wd, an = %wd\n", DIGIT_RADIX(radix), radix->exp, an);
            flint_printf("a = %{ulong*}\n", a, an);
            flint_printf("s  = %{ulong*}\n", s, sn);
            flint_printf("s0 = %{ulong*}\n", s0, sn0);
            flint_abort();
        }

        flint_free(a);
        flint_free(s);
        flint_free(s0);

        radix_clear(radix);
    }

    TEST_FUNCTION_END(state);
}
