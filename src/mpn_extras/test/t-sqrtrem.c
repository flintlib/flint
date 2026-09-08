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

TEST_FUNCTION_START(flint_mpn_sqrtrem, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        mp_ptr a, s, r, s2, r2, t;
        mp_size_t an, sn, rn, rn2;
        int alg, sq, sq2;

        an = 1 + n_randint(state, 60);
        if (n_randint(state, 3) == 0)
            an = 2;
        if (n_randint(state, 10) == 0)
            an = 1 + n_randint(state, 2000);
        if (n_randint(state, 40) == 0)
            an = 2800 + n_randint(state, 1500);
        sn = (an + 1) / 2;

        a = flint_malloc(an * sizeof(mp_limb_t));
        s = flint_malloc((sn + 1) * sizeof(mp_limb_t));
        r = flint_malloc((FLINT_MAX(an, 2) + 1) * sizeof(mp_limb_t));   /* room for an limbs */
        s2 = flint_malloc(sn * sizeof(mp_limb_t));
        r2 = flint_malloc(an * sizeof(mp_limb_t));
        t = flint_malloc(2 * sn * sizeof(mp_limb_t));

        flint_mpn_rrandom(a, state, an);
        if (an == 2 && n_randint(state, 2))
            flint_mpn_urandomb(a, state, 2 * FLINT_BITS - n_randint(state, 64));
        if (a[an - 1] == 0)
            a[an - 1] = 1;
        if (an == 2 && n_randint(state, 20) == 0)
            a[1] = UWORD_MAX;

        /* perfect squares and near squares */
        if (n_randint(state, 3) == 0)
        {
            flint_mpn_rrandom(s2, state, sn);
            flint_mpn_sqr(t, s2, sn);
            if (t[2 * sn - 1] == 0 && an == 2 * sn)
                t[2 * sn - 1] = 1;
            flint_mpn_copyi(a, t, an);
            if (a[an - 1] == 0)
                a[an - 1] = 1;
            if (n_randint(state, 2))
                mpn_add_1(a, a, an, n_randint(state, 3));
        }

        rn2 = mpn_sqrtrem(s2, r2, a, an);

        alg = n_randint(state, 4);
        if (alg == 0)
            rn = flint_mpn_sqrtrem(s, r, a, an);
        else if (alg == 1 && an >= 2)
        {
            _flint_mpn_sqrtrem_newton(s, r, a, an);
            rn = sn + 1;
            while (rn > 0 && r[rn - 1] == 0)
                rn--;
        }
        else
        {
            /* NULL remainder: returns 0 iff perfect square */
            alg = 2;
            rn = flint_mpn_sqrtrem(s, NULL, a, an);
            if (rn != (rn2 != 0))
                TEST_FUNCTION_FAIL("NULL remainder return value: an = %wd, rn = %wd, rn2 = %wd\n", an, rn, rn2);
            rn = rn2;
            flint_mpn_copyi(r, r2, rn2);
            flint_mpn_zero(r + rn2, sn + 1 - rn2);
        }

        if (rn != rn2 || mpn_cmp(s, s2, sn) != 0 || mpn_cmp(r, r2, rn2) != 0
            || !flint_mpn_zero_p(r + rn2, sn + 1 - rn2))
            TEST_FUNCTION_FAIL("alg = %d, an = %wd, rn = %wd, rn2 = %wd\na = %{ulong*}\ns = %{ulong*}\ns2 = %{ulong*}\n",
                alg, an, rn, rn2, a, an, s, sn, s2, sn);

        /* checked exact square root */
        sq = flint_mpn_sqrt(s, a, an);
        sq2 = (rn2 == 0);
        if (sq != sq2 || (sq && mpn_cmp(s, s2, sn) != 0))
            TEST_FUNCTION_FAIL("flint_mpn_sqrt: an = %wd, sq = %d, sq2 = %d\n", an, sq, sq2);

        flint_free(a);
        flint_free(s);
        flint_free(r);
        flint_free(s2);
        flint_free(r2);
        flint_free(t);
    }

    TEST_FUNCTION_END(state);
}
