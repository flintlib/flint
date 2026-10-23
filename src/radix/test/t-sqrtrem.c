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

TEST_FUNCTION_START(radix_sqrtrem, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        nn_ptr a, s, r, s2, r2, P, t;
        slong an, sn;
        ulong one = 1;

        radix_t radix;
        radix_init_randtest(radix, state);

        if (n_randint(state, 100))
            an = 1 + n_randint(state, 30);
        else
            an = 1 + n_randint(state, 300);

        a = flint_malloc(sizeof(ulong) * (an + 2));

        radix_randtest_limbs(a, state, an, radix);
        if (a[an - 1] == 0)
            a[an - 1] = 1 + n_randint(state, LIMB_RADIX(radix) - 1);

        sn = (an + 1) / 2;

        /* Bias towards perfect squares and near-squares to exercise the
           boundary digits. */
        if (n_randint(state, 4) == 0)
        {
            nn_ptr s0 = flint_malloc(sizeof(ulong) * sn);
            radix_randtest_limbs(s0, state, sn, radix);
            if (s0[sn - 1] == 0)
                s0[sn - 1] = 1 + n_randint(state, LIMB_RADIX(radix) - 1);
            radix_mul(a, s0, sn, s0, sn, radix);
            an = 2 * sn;
            MPN_NORM(a, an);
            if (n_randint(state, 2))
            {
                ulong d = 1 + n_randint(state, LIMB_RADIX(radix) - 1);
                a[an] = 0;
                radix_add(a, a, an + 1, &d, 1, radix);
                an = an + 1;
            }
            MPN_NORM(a, an);
            sn = (an + 1) / 2;
            flint_free(s0);
        }

        s = flint_malloc(sizeof(ulong) * sn);
        r = flint_malloc(sizeof(ulong) * (sn + 1));
        s2 = flint_malloc(sizeof(ulong) * sn);
        r2 = flint_malloc(sizeof(ulong) * (sn + 1));
        P = flint_malloc(sizeof(ulong) * (an + 2));
        t = flint_malloc(sizeof(ulong) * (sn + 2));

        radix_sqrtrem(s, r, a, an, radix);

        /* Check s^2 + r == a */
        flint_mpn_zero(P, an + 1);
        radix_mul(P, s, sn, s, sn, radix);
        if (2 * sn < an + 1)
            flint_mpn_zero(P + 2 * sn, an + 1 - 2 * sn);
        radix_add(P, P, an + 1, r, sn + 1, radix);

        if (P[an] != 0 || mpn_cmp(P, a, an) != 0)
        {
            flint_printf("FAIL: s^2 + r != a: radix %wu^%wd, an = %wd\n", DIGIT_RADIX(radix), radix->exp, an);
            flint_printf("a = %{ulong*}\n", a, an);
            flint_printf("s = %{ulong*}\n", s, sn);
            flint_printf("r = %{ulong*}\n", r, sn + 1);
            flint_abort();
        }

        /* Check r < 2s + 1 */
        t[sn] = radix_add(t, s, sn, s, sn, radix);
        radix_add(t, t, sn + 1, &one, 1, radix);

        if (mpn_cmp(r, t, sn + 1) >= 0)
        {
            flint_printf("FAIL: r >= 2s + 1: radix %wu^%wd, an = %wd\n", DIGIT_RADIX(radix), radix->exp, an);
            flint_printf("a = %{ulong*}\n", a, an);
            flint_printf("s = %{ulong*}\n", s, sn);
            flint_printf("r = %{ulong*}\n", r, sn + 1);
            flint_abort();
        }

        /* Cross-check the Newton-Karp-Markstein path against the dispatch
           (which uses mpn conversion for small an), and the r == NULL
           variant. */
        radix_sqrtrem_newton_karp_markstein(s2, r2, a, an, radix);

        if (mpn_cmp(s, s2, sn) != 0 || mpn_cmp(r, r2, sn + 1) != 0)
        {
            flint_printf("FAIL: newton_karp_markstein mismatch: radix %wu^%wd, an = %wd\n", DIGIT_RADIX(radix), radix->exp, an);
            flint_printf("a = %{ulong*}\n", a, an);
            flint_printf("s = %{ulong*}\n", s, sn);
            flint_printf("s2 = %{ulong*}\n", s2, sn);
            flint_abort();
        }

        radix_sqrtrem_newton_karp_markstein(s2, NULL, a, an, radix);

        if (mpn_cmp(s, s2, sn) != 0)
        {
            flint_printf("FAIL: r == NULL mismatch: radix %wu^%wd, an = %wd\n", DIGIT_RADIX(radix), radix->exp, an);
            flint_abort();
        }

        flint_free(a);
        flint_free(s);
        flint_free(r);
        flint_free(s2);
        flint_free(r2);
        flint_free(P);
        flint_free(t);

        radix_clear(radix);
    }

    TEST_FUNCTION_END(state);
}
