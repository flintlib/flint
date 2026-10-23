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
        /* rarely, the Newton square root */
        if (n_randint(state, 1000) == 0)
            an = FLINT_MPN_SQRTREM_NEWTON_CUTOFF + n_randint(state, 300);
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
            else if (n_randint(state, 2) && an > 1)
                mpn_sub_1(a, a, an, 1 + n_randint(state, 2));
            if (a[an - 1] == 0)
                a[an - 1] = 1;
        }

        rn2 = mpn_sqrtrem(s2, r2, a, an);

        alg = n_randint(state, 6);
        if (alg == 0)
            rn = flint_mpn_sqrtrem(s, r, a, an);
        else if (alg == 1 && an >= 2)
        {
            _flint_mpn_sqrtrem_newton(s, r, a, an);
            rn = sn + 1;
            while (rn > 0 && r[rn - 1] == 0)
                rn--;
        }
#if FLINT_BITS == 64
        else if (alg == 3 && an >= 5)
        {
            rn = _flint_mpn_sqrtrem_divconquer(s, r, a, an);
        }
        else if (alg == 4 && an >= 5)
        {
            rn = _flint_mpn_sqrtrem_divconquer(s, NULL, a, an);
            if (rn != (rn2 != 0))
                TEST_FUNCTION_FAIL("dc, NULL remainder return value: an = %wd, rn = %wd, rn2 = %wd\n", an, rn, rn2);
            rn = rn2;
            flint_mpn_copyi(r, r2, rn2);
            flint_mpn_zero(r + rn2, sn + 1 - rn2);
        }
#endif
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

    /* two to four limbs (dedicated code on 64-bit machines): perfect and
       near squares, remainders close to 2s, normalization boundaries and
       roots close to a limb boundary */
    for (iter = 0; iter < 100000 * flint_test_multiplier(); iter++)
    {
        mp_limb_t a[4], s[2], r[5], s2[2], r2[4], t[4], u[2], v[3];
        mp_size_t an, sn, rn, rn2;
        int i, kind;

        an = 2 + n_randint(state, 3);
        sn = (an + 1) / 2;
        kind = n_randint(state, 8);

        for (i = 0; i < an; i++)
            a[i] = n_randlimb(state);

        if (kind == 1)
        {
            a[an - 1] >>= n_randint(state, FLINT_BITS);
        }
        else if (kind == 2)
        {
            for (i = 0; i < an; i++)
                a[i] = n_randtest(state);
        }
        else if (kind == 3)
        {
            /* leading limbs all ones or nearly */
            for (i = an - 1; i >= 0 && n_randint(state, 4) != 0; i--)
                a[i] = UWORD_MAX - n_randint(state, 3);
        }
        else if (kind == 4 || kind == 5)
        {
            /* u^2 + k or (u + 1)^2 - 1 - k = u^2 + 2u - k for small k,
               with u of sn limbs chosen so that u^2 has an limbs */
            u[0] = n_randlimb(state);
            u[1] = (sn == 2) ? n_randlimb(state) : 0;
            if (an == 2)
                u[0] |= UWORD(1) << (FLINT_BITS / 2 + n_randint(state, FLINT_BITS / 2));
            else if (an == 3)
                u[1] = (u[1] >> (FLINT_BITS / 2)) | 1;
            else
                u[1] |= UWORD(1) << (FLINT_BITS / 2 + n_randint(state, FLINT_BITS / 2));
            if (n_randint(state, 4) == 0)
            {
                /* root close to a power of two or to the limb boundary */
                u[sn - 1] = (an == 3) ? (UWORD(1) << (FLINT_BITS / 2)) - 1 : UWORD_MAX;
                u[0] = (sn == 2) ? UWORD_MAX - n_randint(state, 3) : u[0];
            }
            flint_mpn_sqr(t, u, sn);
            flint_mpn_copyi(a, t, an);
            if (kind == 4)
            {
                mpn_add_1(a, a, an, n_randint(state, 3));
            }
            else
            {
                v[sn] = mpn_lshift(v, u, sn, 1);
                mpn_add(a, a, an, v, sn + 1);
                mpn_sub_1(a, a, an, n_randint(state, 3));
            }
        }
        else if (kind == 6)
        {
            /* top limb near a power of two (normalization shift boundaries) */
            a[an - 1] = UWORD(1) << n_randint(state, FLINT_BITS);
            if (n_randint(state, 2))
                a[an - 1] -= (a[an - 1] > 1);
            if (n_randint(state, 2))
                a[an - 2] = n_randint(state, 2) ? 0 : UWORD_MAX;
        }
        else if (kind == 7)
        {
            /* top two limbs a square or one below */
            u[0] = n_randlimb(state) >> n_randint(state, FLINT_BITS / 2);
            if (n_randint(state, 4) == 0)
                u[0] = UWORD_MAX - n_randint(state, 2);
            umul_ppmm(a[an - 1], a[an - 2], u[0], u[0]);
            if (n_randint(state, 2))
                sub_ddmmss(a[an - 1], a[an - 2], a[an - 1], a[an - 2], 0, 1);
        }

        if (a[an - 1] == 0)
            a[an - 1] = 1 + n_randint(state, 2);

        rn2 = mpn_sqrtrem(s2, r2, a, an);
        for (i = rn2; i < sn + 1; i++)
            r2[i] = 0;

        for (i = 0; i < 5; i++)
            r[i] = n_randlimb(state);
        rn = flint_mpn_sqrtrem(s, r, a, an);

        if (rn != rn2 || mpn_cmp(s, s2, sn) != 0 || mpn_cmp(r, r2, sn + 1) != 0)
            TEST_FUNCTION_FAIL("small: an = %wd, kind = %d, rn = %wd, rn2 = %wd\na = %{ulong*}\ns = %{ulong*}\ns2 = %{ulong*}\nr = %{ulong*}\nr2 = %{ulong*}\n",
                an, kind, rn, rn2, a, an, s, sn, s2, sn, r, sn + 1, r2, sn + 1);

        s[0] = s[sn - 1] = 0;
        rn = flint_mpn_sqrtrem(s, NULL, a, an);
        if (rn != (rn2 != 0) || mpn_cmp(s, s2, sn) != 0)
            TEST_FUNCTION_FAIL("small, NULL remainder: an = %wd, kind = %d, rn = %wd, rn2 = %wd\na = %{ulong*}\n",
                an, kind, rn, rn2, a, an);
    }

    TEST_FUNCTION_END(state);
}
