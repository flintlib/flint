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

/* n_ll_barrett_precomp_preinv, n_ll_mod_barrett, n_lll_mod_barrett */

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

TEST_FUNCTION_START(n_ll_mod_barrett, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * flint_test_multiplier(); iter++)
    {
        ulong n, ninv, qhi, qlo, b1, b2, a[3], q[3], r, s;
        int bits = 2 + n_randint(state, FLINT_BITS - 2);

        n = n_randbits(state, bits);    /* 2 <= n < 2^(FLINT_BITS - 1) */
        if (n_randint(state, 8) == 0)
            n = UWORD(1) << (bits - 1);
        if (n_randint(state, 8) == 0)
            n = (UWORD(1) << (FLINT_BITS - 1)) - 1 - n_randint(state, 4);
        if (n < 2)
            n = 2;
        /* full-word modulus: only the constants, which the CRT leaves also
           use for such moduli; the reductions require n < 2^(FLINT_BITS - 1) */
        if (n_randint(state, 8) == 0)
            n = n_randtest(state) | (UWORD(1) << (FLINT_BITS - 1));

        ninv = n_preinvert_limb(n);
        n_ll_barrett_precomp_preinv(&qhi, &qlo, &b1, &b2, n, ninv);

        /* the constants */
        a[0] = 0; a[1] = 0; a[2] = 1;
        mpn_divrem_1(q, 0, a, 3, n);
        if (qhi != q[1] || qlo != q[0] || q[2] != 0)
            TEST_FUNCTION_FAIL("floor(B^2 / n)\nn = %wu, got (%wu, %wu), expected (%wu, %wu)\n",
                               n, qhi, qlo, q[1], q[0]);
        a[0] = 0; a[1] = 1;
        r = mpn_mod_1(a, 2, n);
        a[0] = 0; a[1] = 0; a[2] = 1;
        s = mpn_mod_1(a, 3, n);
        if (b1 != r || b2 != s)
            TEST_FUNCTION_FAIL("B mod n, B^2 mod n\nn = %wu, got (%wu, %wu), expected (%wu, %wu)\n",
                               n, b1, b2, r, s);

        if (n >= (UWORD(1) << (FLINT_BITS - 1)))
            continue;

        /* two words */
        a[0] = randtest_word(state);
        a[1] = randtest_word(state);
        r = n_ll_mod_barrett(a[1], a[0], n, qhi, qlo);
        s = mpn_mod_1(a, 2, n);
        if (r != s)
            TEST_FUNCTION_FAIL("n_ll_mod_barrett\nn = %wu, a = (%wu, %wu), got %wu, expected %wu\n",
                               n, a[1], a[0], r, s);

        /* three words */
        a[0] = randtest_word(state);
        a[1] = randtest_word(state);
        a[2] = randtest_word(state);
        r = n_lll_mod_barrett(a[2], a[1], a[0], n, qhi, qlo, b1, b2);
        s = mpn_mod_1(a, 3, n);
        if (r != s)
            TEST_FUNCTION_FAIL("n_lll_mod_barrett\nn = %wu, a = (%wu, %wu, %wu), got %wu, expected %wu\n",
                               n, a[2], a[1], a[0], r, s);
    }

    TEST_FUNCTION_END(state);
}
