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

TEST_FUNCTION_START(radix_divrem, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        radix_t radix;
        nn_ptr a, b, q, q2, r, bq, bqr;
        slong an, bn, qn, rn;
        int algorithm, aliasing;

        radix_init_randtest(radix, state);

        slong N = n_randint(state, 100) ? 6 : 100;
        bn = 1 + n_randint(state, N);

        an = bn + n_randint(state, N);
        qn = an - bn + 1;
        rn = bn;

        a = flint_malloc(an * sizeof(ulong));
        b = flint_malloc(bn * sizeof(ulong));
        q = flint_malloc(qn * sizeof(ulong));
        q2 = flint_malloc(qn * sizeof(ulong));
        r = flint_malloc(rn * sizeof(ulong));
        bq = flint_malloc(an * sizeof(ulong));
        bqr = flint_malloc(an * sizeof(ulong));

        radix_randtest_limbs(a, state, an, radix);
        radix_randtest_limbs(b, state, bn - 1, radix);
        b[bn - 1] = n_randint(state, LIMB_RADIX(radix));
        b[bn - 1] = FLINT_MAX(b[bn - 1], 1);

        radix_randtest_limbs(q, state, qn, radix);
        radix_randtest_limbs(r, state, rn, radix);

        algorithm = n_randint(state, 3);

        void (*divrem_func)(nn_ptr q, nn_ptr r,
                           nn_srcptr a, slong an,
                           nn_srcptr b, slong bn,
                           const radix_t radix);

        if (algorithm == 0)
            divrem_func = radix_divrem_via_mpn;
        else if (algorithm == 1)
            divrem_func = radix_divrem_newton;
        else
            divrem_func = radix_divrem;

        aliasing = n_randint(state, 7);

        if (aliasing == 0)
        {
            divrem_func(q, r, a, an, b, bn, radix);
        }
        else if (aliasing == 1)
        {
            q = flint_realloc(q, an * sizeof(ulong));
            flint_mpn_copyi(q, a, an);
            divrem_func(q, r, q, an, b, bn, radix);
        }
        else if (aliasing == 2)
        {
            r = flint_realloc(r, an * sizeof(ulong));
            flint_mpn_copyi(r, a, an);
            divrem_func(q, r, r, an, b, bn, radix);
        }
        else if (aliasing == 3)
        {
            q = flint_realloc(q, FLINT_MAX(qn, bn) * sizeof(ulong));
            flint_mpn_copyi(q, b, bn);
            divrem_func(q, r, a, an, q, bn, radix);
        }
        else if (aliasing == 4)
        {
            flint_mpn_copyi(r, b, bn);
            divrem_func(q, r, a, an, r, bn, radix);
        }
        else if (aliasing == 5)
        {
            q = flint_realloc(q, an * sizeof(ulong));
            flint_mpn_copyi(q, a, an);
            flint_mpn_copyi(r, b, bn);
            divrem_func(q, r, q, an, r, bn, radix);
        }
        else if (aliasing == 6)
        {
            q = flint_realloc(q, FLINT_MAX(qn, bn) * sizeof(ulong));
            flint_mpn_copyi(q, b, bn);
            r = flint_realloc(r, an * sizeof(ulong));
            flint_mpn_copyi(r, a, an);
            divrem_func(q, r, r, an, q, bn, radix);
        }

        radix_mulmid(bq, q, qn, b, bn, 0, an, radix);
        radix_add(bqr, bq, an, r, rn, radix);

        if (mpn_cmp(bqr, a, an) != 0 || mpn_cmp(r, b, rn) > 0)
        {
            flint_printf("FAIL: radix_divrem\n");
            flint_printf("algorithm = %d\n", algorithm);
            flint_printf("aliasing = %d\n", aliasing);
            flint_printf("a = %{ulong*}\n", a, an);
            flint_printf("b = %{ulong*}\n", b, bn);
            flint_printf("q = %{ulong*}\n", q, qn);
            flint_printf("r = %{ulong*}\n", r, rn);
            flint_printf("bq = %{ulong*}\n", bq, an);
            flint_printf("bqr = %{ulong*}\n", bqr, an);
            flint_abort();
        }

        aliasing = n_randint(state, 2);

        if (aliasing == 0)
        {
            divrem_func(q2, NULL, a, an, b, bn, radix);
        }
        else
        {
            q2 = flint_realloc(q2, an * sizeof(ulong));
            flint_mpn_copyi(q2, a, an);
            divrem_func(q2, NULL, q2, an, b, bn, radix);
        }

        if (mpn_cmp(q, q2, qn) != 0)
        {
            flint_printf("FAIL: radix_divrem (r == NULL)\n");
            flint_printf("algorithm = %d\n", algorithm);
            flint_printf("aliasing = %d\n", aliasing);
            flint_printf("a = %{ulong*}\n", a, an);
            flint_printf("b = %{ulong*}\n", b, bn);
            flint_printf("q = %{ulong*}\n", q, qn);
            flint_printf("q2 = %{ulong*}\n", q2, qn);
            flint_abort();
        }

        radix_clear(radix);

        flint_free(a);
        flint_free(b);
        flint_free(q);
        flint_free(q2);
        flint_free(r);
        flint_free(bq);
        flint_free(bqr);
    }

    TEST_FUNCTION_END(state);
}
