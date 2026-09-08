/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"

TEST_FUNCTION_START(flint_mpn_divexact, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        mp_ptr a, b, q, q2;
        mp_size_t an, bn, n;
        int alg, exact;

        bn = 1 + n_randint(state, 30);
        n = 1 + n_randint(state, 40);
        if (n_randint(state, 20) == 0)
        {
            bn = 1 + n_randint(state, 800);
            n = 1 + n_randint(state, 800);
        }

        b = flint_malloc(bn * sizeof(mp_limb_t));
        q = flint_malloc((n + 1) * sizeof(mp_limb_t));
        q2 = flint_malloc(n * sizeof(mp_limb_t));
        a = flint_malloc((n + bn) * sizeof(mp_limb_t));

        flint_mpn_rrandom(b, state, bn);
        flint_mpn_rrandom(q2, state, n);
        if (b[bn - 1] == 0)
            b[bn - 1] = 1;
        if (n_randint(state, 4) == 0)
            b[0] <<= n_randint(state, FLINT_BITS);
        if (n_randint(state, 8) == 0)
            flint_mpn_zero(b, n_randint(state, bn));
        if (b[bn - 1] == 0)
            b[bn - 1] = 1;
        if (q2[n - 1] == 0)
            q2[n - 1] = 1;

        /* a = q2 b with an = n + bn - 1 or n + bn limbs */
        if (n >= bn)
            flint_mpn_mul(a, q2, n, b, bn);
        else
            flint_mpn_mul(a, b, bn, q2, n);
        an = n + bn - (a[n + bn - 1] == 0);

        alg = n_randint(state, 4);
        if (alg == 0)
            flint_mpn_divexact(q, a, an, b, bn);
        else if (alg == 1)
            _flint_mpn_divexact_hensel(q, a, an, b, bn);
        else if (alg == 3)
        {
            flint_mpn_divexact_preinv_t pre;
            flint_mpn_divexact_preinv_init(pre, b, bn);
            flint_mpn_divexact_preinv(q, a, an, pre);
            flint_mpn_divexact_preinv_clear(pre);
        }
        else
        {
            exact = flint_mpn_div(q, a, an, b, bn);
            if (!exact)
                TEST_FUNCTION_FAIL("flint_mpn_div returned 0 for exact division: an = %wd, bn = %wd\n", an, bn);
        }

        /* the quotient has an - bn + 1 limbs, possibly one more than n */
        if (mpn_cmp(q, q2, n) != 0 || (an - bn + 1 > n && q[n] != 0))
            TEST_FUNCTION_FAIL("alg = %d, an = %wd, bn = %wd, n = %wd\na = %{ulong*}\nb = %{ulong*}\nq = %{ulong*}\nq2 = %{ulong*}\n",
                alg, an, bn, n, a, an, b, bn, q, n, q2, n);

        /* non-exact test for flint_mpn_div */
        if (alg == 2 && bn > 1)
        {
            mp_limb_t cy = mpn_add_1(a, a, an, 1 + n_randint(state, 5));
            if (cy && an < n + bn)
            {
                a[an] = cy;
                an++;
                cy = 0;
            }
            exact = cy ? 0 : flint_mpn_div(q, a, an, b, bn);
            if (exact)
                TEST_FUNCTION_FAIL("flint_mpn_div returned 1 for inexact division: an = %wd, bn = %wd\n", an, bn);
        }

        flint_free(a);
        flint_free(b);
        flint_free(q);
        flint_free(q2);
    }

    TEST_FUNCTION_END(state);
}
