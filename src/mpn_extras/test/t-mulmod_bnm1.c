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

/* canonical reduction mod B^rn - 1 for comparison */
static void
canon_bnm1(nn_ptr r, nn_srcptr v, mp_size_t vn, mp_size_t rn)
{
    mp_size_t c, i;
    ulong cy;

    c = FLINT_MIN(vn, rn);
    flint_mpn_copyi(r, v, c);
    flint_mpn_zero(r + c, rn - c);
    v += c; vn -= c;
    while (vn > 0)
    {
        c = FLINT_MIN(vn, rn);
        cy = mpn_add(r, r, rn, v, c);
        while (cy)
            cy = mpn_add_1(r, r, rn, cy);
        v += c; vn -= c;
    }
    for (i = 0; i < rn; i++)
        if (r[i] != UWORD_MAX)
            return;
    flint_mpn_zero(r, rn);
}

TEST_FUNCTION_START(flint_mpn_mulmod_bnm1, state)
{
    for (slong iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        mp_size_t n = 1 + n_randint(state, n_randint(state, 20) == 0 ?
                                           900 : 120);
        mp_size_t rn = flint_mpn_mulmod_bnm1_next_size(n);
        mp_size_t an = 1 + n_randint(state, rn);
        mp_size_t bn = 1 + n_randint(state, rn);
        int sqr = n_randint(state, 4) == 0;
        nn_ptr a, b, r, ref, v, tp;
        mp_size_t i;

        if (rn < n)
            TEST_FUNCTION_FAIL("next_size: n = %wd, rn = %wd\n", n, rn);

        a = FLINT_ARRAY_ALLOC(rn, ulong);
        b = FLINT_ARRAY_ALLOC(rn, ulong);
        r = FLINT_ARRAY_ALLOC(2 * rn, ulong);
        ref = FLINT_ARRAY_ALLOC(rn, ulong);
        v = FLINT_ARRAY_ALLOC(2 * rn, ulong);
        tp = FLINT_ARRAY_ALLOC(flint_mpn_mulmod_bnm1_itch(rn), ulong);

        flint_mpn_rrandom(a, state, an);
        if (n_randint(state, 20) == 0)
            for (i = 0; i < an; i++)
                a[i] = UWORD_MAX;
        if (sqr)
        {
            bn = an;
            flint_mpn_copyi(b, a, bn);
        }
        else
            flint_mpn_rrandom(b, state, bn);

        flint_mpn_mulmod_bnm1(r, rn, a, an, sqr ? a : b, bn, tp);

        /* representative < B^rn congruent to the true product */
        if (an >= bn)
            mpn_mul(v, a, an, sqr ? a : b, bn);
        else
            mpn_mul(v, sqr ? a : b, bn, a, an);
        canon_bnm1(ref, v, an + bn, rn);
        canon_bnm1(v, r, rn, rn);

        if (mpn_cmp(v, ref, rn) != 0)
            TEST_FUNCTION_FAIL("n = %wd, rn = %wd, an = %wd, bn = %wd, "
                               "sqr = %d\n", n, rn, an, bn, sqr);

        flint_free(a); flint_free(b); flint_free(r);
        flint_free(ref); flint_free(v); flint_free(tp);
    }

    TEST_FUNCTION_END(state);
}
