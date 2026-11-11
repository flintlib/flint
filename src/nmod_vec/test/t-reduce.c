/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_vec.h"

/* a word to reduce, often at the corners of the SIMD kernels: q n + d with
   q at the top and d small, low half all ones, a power of 2, ... */
static ulong
_red_word(flint_rand_t state, ulong n)
{
    ulong t, d;

    switch (n_randint(state, 7))
    {
        case 0:
            return UWORD_MAX - n_randint(state, 3);
        case 1:
            t = (UWORD_MAX / n - n_randint(state, 3)) * n;
            d = n_randint(state, 3);
            return (t + d >= t) ? t + d : t;
        case 2:
            return n_randlimb(state) | ((UWORD(1) << (FLINT_BITS / 2)) - 1);
        case 3:
            return UWORD(1) << n_randint(state, FLINT_BITS);
        case 4:
            return n_randint(state, 3);
        case 5:
            return n_randtest(state);
        default:
            return n_randlimb(state);
    }
}

/* the moduli 2^(nu - 1), 2^(nu - 1) + 1, 2^nu - 1 and a random one */
static ulong
_red_modulus(flint_rand_t state, int nu, int k)
{
    ulong lo = UWORD(1) << (nu - 1);
    ulong hi = (nu == FLINT_BITS) ? UWORD_MAX : (UWORD(1) << nu) - 1;

    if (k == 0 || lo == hi)
        return lo;
    else if (k == 1)
        return lo + 1;
    else if (k == 2)
        return hi;
    else
        return lo + n_randint(state, hi - lo + 1);
}

/* checks vec[i] = a[i] mod n for i < len (lazy: congruent and below 2n),
   and that the guard words of r around [off, off + len) are unchanged */
static int
_red_check(nn_srcptr r, nn_srcptr a, slong off, slong len, ulong n, ulong guard,
           int lazy)
{
    slong i;

    for (i = 0; i < off; i++)
        if (r[i] != guard)
            return 0;
    for (i = 0; i < len; i++)
        if (lazy ? (r[off + i] / 2 >= n || r[off + i] % n != a[i] % n)
                 : r[off + i] != a[i] % n)
            return 0;
    for (i = off + len; i < off + len + 8; i++)
        if (r[i] != guard)
            return 0;
    return 1;
}

#define RED_MAXLEN 256

TEST_FUNCTION_START(nmod_vec_reduce, state)
{
    const ulong guard = UWORD_MAX / 3;
    nn_ptr a, r, a0;
    slong iter, i, len, off;
    int nu, k, inplace, lazy;

    /* 64-byte aligned buffers with room for offsets and guard words */
    a = flint_aligned_alloc(64, (RED_MAXLEN + 16) * sizeof(ulong));
    r = flint_aligned_alloc(64, (RED_MAXLEN + 16) * sizeof(ulong));
    a0 = flint_malloc(RED_MAXLEN * sizeof(ulong));

    /* every bit size, lengths 0 to 3 vectors of 8 words plus a tail, and
       around the peeling threshold of 128; offsets 0, 1, 2 and 7 words
       from a 64-byte boundary; in place and not; _nmod_vec_reduce_lazy up
       to 62 bits */
    for (nu = 1; nu <= FLINT_BITS; nu++)
    for (k = 0; k < 4; k++)
    {
        nmod_t mod;
        ulong n = _red_modulus(state, nu, k);
        nmod_init(&mod, n);

        for (len = 0; len <= 140; len += (len < 35 || len >= 120) ? 1 : 17)
        for (off = 0; off < 8; off += (off < 2) ? 1 : 5)
        for (inplace = 0; inplace < 2; inplace++)
        for (lazy = 0; lazy < 1 + (nu <= 62); lazy++)
        {
            nn_ptr out = inplace ? a : r;

            for (i = 0; i < RED_MAXLEN + 16; i++)
                a[i] = r[i] = guard;
            for (i = 0; i < len; i++)
                a0[i] = a[off + i] = _red_word(state, n);

            if (lazy)
                _nmod_vec_reduce_lazy(out + off, a + off, len, mod);
            else
                _nmod_vec_reduce(out + off, a + off, len, mod);
            if (!_red_check(out, a0, off, len, n, guard, lazy))
                TEST_FUNCTION_FAIL("nu = %d, n = %wu, len = %wd, off = %wd, in place %d, lazy %d\n",
                                   nu, n, len, off, inplace, lazy);
        }
    }

    /* random moduli and lengths, the input and the output at different
       offsets */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        nmod_t mod;
        ulong n = n_randtest_not_zero(state);
        slong offa = n_randint(state, 8), offr = n_randint(state, 8);

        len = n_randint(state, RED_MAXLEN);
        nmod_init(&mod, n);

        for (i = 0; i < RED_MAXLEN + 16; i++)
            a[i] = r[i] = guard;
        for (i = 0; i < len; i++)
            a0[i] = a[offa + i] = _red_word(state, n);

        lazy = NMOD_BITS(mod) <= 62 && n_randint(state, 2);
        if (lazy)
            _nmod_vec_reduce_lazy(r + offr, a + offa, len, mod);
        else
            _nmod_vec_reduce(r + offr, a + offa, len, mod);
        if (!_red_check(r, a0, offr, len, n, guard, lazy))
            TEST_FUNCTION_FAIL("n = %wu, len = %wd, offsets %wd, %wd\n", n, len, offa, offr);
    }

    flint_aligned_free(a);
    flint_aligned_free(r);
    flint_free(a0);

    TEST_FUNCTION_END(state);
}
