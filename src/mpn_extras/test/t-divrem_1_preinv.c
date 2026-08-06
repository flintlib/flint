/*
    Copyright (C) 2013 William Hart
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "mpn_extras.h"

TEST_FUNCTION_START(flint_mpn_divrem_1_preinv, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        mp_ptr a, q1, q2;
        ulong d, dinv, norm, r1, r2;
        int alias = n_randint(state, 2);
        slong n = 1 + n_randint(state, 10);

        a = flint_malloc(sizeof(mp_limb_t) * n);
        q1 = flint_malloc(sizeof(mp_limb_t) * n);
        q2 = flint_malloc(sizeof(mp_limb_t) * n);

        d = n_randtest_not_zero(state);
        norm = flint_clz(d);
        dinv = n_preinvert_limb_prenorm(d << norm);

        flint_mpn_rrandom(a, state, n);

        if (alias)
        {
            flint_mpn_copyi(q1, a, n);
            r1 = flint_mpn_divrem_1_preinv(q1, q1, n, d, dinv, norm);
        }
        else
        {
            r1 = flint_mpn_divrem_1_preinv(q1, a, n, d, dinv, norm);
        }
        r2 = mpn_divrem_1(q2, 0, a, n, d);

        if (mpn_cmp(q1, q2, n) || r1 != r2)
        {
            flint_printf("FAIL: flint_mpn_divrem_1_preinv\n");
            flint_printf("d = %wu, n = %wd, alias = %d\n", d, n, alias);
            flint_printf("a = %{ulong*}\n", a, n);
            flint_printf("q1 = %{ulong*}\n", q1, n);
            flint_printf("r1 = %wu\n", r1);
            flint_printf("q2 = %{ulong*}\n", q2, n);
            flint_printf("r2 = %wu\n", r2);
            flint_abort();
        }

        flint_free(a);
        flint_free(q1);
        flint_free(q2);
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        mp_ptr a, q1, q2;
        ulong d, dinv, norm, r1, r2;
        int alias = n_randint(state, 2);

        a = flint_malloc(sizeof(mp_limb_t) * 2);
        q1 = flint_malloc(sizeof(mp_limb_t) * 2);
        q2 = flint_malloc(sizeof(mp_limb_t) * 2);

        d = n_randtest_not_zero(state);
        norm = flint_clz(d);
        dinv = n_preinvert_limb_prenorm(d << norm);

        flint_mpn_rrandom(a, state, 2);

        if (alias)
        {
            flint_mpn_copyi(q1, a, 2);

            if (norm == 0)
                r1 = flint_mpn_divrem_2_1_preinv_norm(q1, q1, d, dinv);
            else
                r1 = flint_mpn_divrem_2_1_preinv_unnorm(q1, q1, d, dinv, norm);
        }
        else
        {
            if (norm == 0)
                r1 = flint_mpn_divrem_2_1_preinv_norm(q1, a, d, dinv);
            else
                r1 = flint_mpn_divrem_2_1_preinv_unnorm(q1, a, d, dinv, norm);
        }
        r2 = mpn_divrem_1(q2, 0, a, 2, d);

        if (mpn_cmp(q1, q2, 2) || r1 != r2)
        {
            flint_printf("FAIL: flint_mpn_divrem_2_1_preinv_norm/unnorm\n");
            flint_printf("d = %wu, alias = %d\n", d, alias);
            flint_printf("a = %{ulong*}\n", a, 2);
            flint_printf("q1 = %{ulong*}\n", q1, 2);
            flint_printf("r1 = %wu\n", r1);
            flint_printf("q2 = %{ulong*}\n", q2, 2);
            flint_printf("r2 = %wu\n", r2);
            flint_abort();
        }

        flint_free(a);
        flint_free(q1);
        flint_free(q2);
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        mp_ptr a, q1, q2;
        ulong d, dinv, norm, r1, r2;
        int alias = n_randint(state, 2);

        a = flint_malloc(sizeof(mp_limb_t) * 3);
        q1 = flint_malloc(sizeof(mp_limb_t) * 3);
        q2 = flint_malloc(sizeof(mp_limb_t) * 3);

        d = n_randtest_not_zero(state);
        norm = flint_clz(d);
        dinv = n_preinvert_limb_prenorm(d << norm);

        flint_mpn_rrandom(a, state, 3);

        if (alias)
        {
            flint_mpn_copyi(q1, a, 3);
            if (norm == 0)
                r1 = flint_mpn_divrem_3_1_preinv_norm(q1, q1, d, dinv);
            else
                r1 = flint_mpn_divrem_3_1_preinv_unnorm(q1, q1, d, dinv, norm);
        }
        else
        {
            if (norm == 0)
                r1 = flint_mpn_divrem_3_1_preinv_norm(q1, a, d, dinv);
            else
                r1 = flint_mpn_divrem_3_1_preinv_unnorm(q1, a, d, dinv, norm);
        }
        r2 = mpn_divrem_1(q2, 0, a, 3, d);

        if (mpn_cmp(q1, q2, 3) || r1 != r2)
        {
            flint_printf("FAIL: flint_mpn_divrem_3_1_preinv_norm/unnorm\n");
            flint_printf("d = %wu, alias = %d\n", d, alias);
            flint_printf("a = %{ulong*}\n", a, 3);
            flint_printf("q1 = %{ulong*}\n", q1, 3);
            flint_printf("r1 = %wu\n", r1);
            flint_printf("q2 = %{ulong*}\n", q2, 3);
            flint_printf("r2 = %wu\n", r2);
            flint_abort();
        }

        flint_free(a);
        flint_free(q1);
        flint_free(q2);
    }

    TEST_FUNCTION_END(state);
}
