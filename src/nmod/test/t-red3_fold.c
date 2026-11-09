/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "mpn_extras.h"
#include "nmod.h"

static ulong
randtest_word(flint_rand_t state)
{
    switch (n_randint(state, 4))
    {
        case 0: return UWORD_MAX - n_randint(state, 4);
        case 1: return n_randint(state, 4);
        default: return n_randtest(state);
    }
}

/* a modulus n >= 2 of nu bits; often full-word, close to 2^(FLINT_BITS - 1)
   or to 2^FLINT_BITS */
static ulong
randtest_modulus(flint_rand_t state)
{
    ulong n;

    switch (n_randint(state, 8))
    {
        case 0:
            n = (UWORD(1) << (FLINT_BITS - 1)) + n_randint(state, UWORD(1) << n_randint(state, FLINT_BITS - 1));
            break;
        case 1:
            n = UWORD_MAX - n_randint(state, UWORD(1) << n_randint(state, FLINT_BITS - 1));
            break;
        case 2:
        case 3:
            n = n_randtest(state) | (UWORD(1) << (FLINT_BITS - 1));
            break;
        default:
            n = n_randbits(state, 2 + n_randint(state, FLINT_BITS - 1));
    }

    return FLINT_MAX(n, 2);
}

/* NMOD_RED3_FOLD_FULLWORD, NMOD_RED3_FOLD_FULLWORD_NOSUB with
   _nmod_red3_fold_can_skip_sub, and _nmod_red3_precomp, against mpn_mod_1 */
TEST_FUNCTION_START(nmod_red3_fold, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * flint_test_multiplier(); iter++)
    {
        nmod_t mod;
        nmod_red3_precomp_struct P;
        ulong n, ymax, a[3], r, s, h, l, t1, t0;
        int ok;

        n = randtest_modulus(state);
        nmod_init(&mod, n);

        /* a bound on the top word: small, random or none */
        switch (n_randint(state, 4))
        {
            case 0: ymax = n_randint(state, 64); break;
            case 1: ymax = n_randtest(state) >> n_randint(state, FLINT_BITS); break;
            case 2: ymax = UWORD_MAX; break;
            default: ymax = n_randint(state, UWORD(1) << n_randint(state, FLINT_BITS / 2 + 4));
        }

        _nmod_red3_precomp_init(&P, ymax, mod);

        a[0] = randtest_word(state);
        a[1] = randtest_word(state);
        a[2] = (ymax == UWORD_MAX) ? randtest_word(state) : ymax - n_randint(state, FLINT_MIN(ymax, 4) + 1);
        if (n_randint(state, 4) == 0)
            a[0] = a[1] = UWORD_MAX;
        s = mpn_mod_1(a, 3, n);

        r = _nmod_red3_precomp(a[2], a[1], a[0], &P, mod);
        if (r != s)
            TEST_FUNCTION_FAIL("_nmod_red3_precomp\nn = %wu, ymax = %wu, a = (%wu, %wu, %wu), got %wu, expected %wu\n",
                               n, ymax, a[2], a[1], a[0], r, s);

        if (mod.norm != 0)
            continue;

        /* any three words */
        a[2] = randtest_word(state);
        s = mpn_mod_1(a, 3, n);
        NMOD_RED3_FOLD_FULLWORD(r, a[2], a[1], a[0], P.b1, P.b2, mod);
        if (r != s)
            TEST_FUNCTION_FAIL("NMOD_RED3_FOLD_FULLWORD\nn = %wu, a = (%wu, %wu, %wu), got %wu, expected %wu\n",
                               n, a[2], a[1], a[0], r, s);

        /* the predicate is exact: the fold of (ymax, B - 1, B - 1) has a high
           word below n if and only if it holds */
        if (ymax != UWORD_MAX)
        {
            ulong c;

            ok = _nmod_red3_fold_can_skip_sub(ymax, P.b1, P.b2, mod);
            umul_ppmm(h, l, ymax, P.b2);
            umul_ppmm(t1, t0, UWORD_MAX, P.b1);
            add_sssaaaaaa(c, h, l, UWORD(0), h, l, UWORD(0), t1, t0);
            add_sssaaaaaa(c, h, l, c, h, l, UWORD(0), UWORD(0), UWORD_MAX);
            if (ok != (c == 0 && h < n))
                TEST_FUNCTION_FAIL("_nmod_red3_fold_can_skip_sub\nn = %wu, ymax = %wu, got %d\n",
                                   n, ymax, ok);

            if (ok)
            {
                a[2] = n_randint(state, ymax + 1);
                if (n_randint(state, 2))
                    a[2] = ymax;
                s = mpn_mod_1(a, 3, n);
                NMOD_RED3_FOLD_FULLWORD_NOSUB(r, a[2], a[1], a[0], P.b1, P.b2, mod);
                if (r != s)
                    TEST_FUNCTION_FAIL("NMOD_RED3_FOLD_FULLWORD_NOSUB\nn = %wu, ymax = %wu, a = (%wu, %wu, %wu), got %wu, expected %wu\n",
                                       n, ymax, a[2], a[1], a[0], r, s);
            }
        }
    }

    TEST_FUNCTION_END(state);
}
