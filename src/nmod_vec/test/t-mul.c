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
#include "nmod.h"
#include "nmod_vec.h"

/* a residue, often at the corners: n - 1 - d for small d, small, a power
   of 2 reduced, or random */
static ulong
_mul_residue(flint_rand_t state, ulong n)
{
    ulong d;

    switch (n_randint(state, 5))
    {
        case 0:
            d = n_randint(state, 3);
            return (d < n) ? n - 1 - d : n - 1;
        case 1:
            return n_randint(state, 3) % n;
        case 2:
            return (UWORD(1) << n_randint(state, FLINT_BITS)) % n;
        default:
            return n_randint(state, n);
    }
}

/* the moduli 2^(nu - 1), 2^(nu - 1) + 1, 2^nu - 1 and a random one */
static ulong
_mul_modulus(flint_rand_t state, int nu, int k)
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

/* r[off + i] = a[i] b[i] mod n for i < len, guard words around */
static int
_mul_check(nn_srcptr r, nn_srcptr a, nn_srcptr b, slong off, slong len, nmod_t mod, ulong guard)
{
    slong i;

    for (i = 0; i < off; i++)
        if (r[i] != guard)
            return 0;
    for (i = 0; i < len; i++)
        if (r[off + i] != nmod_mul(a[i], b[i], mod))
            return 0;
    for (i = off + len; i < off + len + 8; i++)
        if (r[i] != guard)
            return 0;
    return 1;
}

#define MUL_MAXLEN 256

TEST_FUNCTION_START(nmod_vec_mul, state)
{
    const ulong guard = UWORD_MAX / 3;
    nn_ptr a, b, r, a0, b0;
    slong iter, i, len, off;
    int nu, k, alias;

    a = flint_aligned_alloc(64, (MUL_MAXLEN + 16) * sizeof(ulong));
    b = flint_aligned_alloc(64, (MUL_MAXLEN + 16) * sizeof(ulong));
    r = flint_aligned_alloc(64, (MUL_MAXLEN + 16) * sizeof(ulong));
    a0 = flint_malloc(MUL_MAXLEN * sizeof(ulong));
    b0 = flint_malloc(MUL_MAXLEN * sizeof(ulong));

    /* every bit size, lengths 0 to 3 vectors of 8 words plus a tail, and
       around the peeling threshold of 128; offsets 0, 1, 2 and 7 words
       from a 64-byte boundary; the output apart, equal to the first input
       or to the second */
    for (nu = 1; nu <= FLINT_BITS; nu++)
    for (k = 0; k < 4; k++)
    {
        nmod_t mod;
        ulong n = _mul_modulus(state, nu, k);
        nmod_init(&mod, n);

        for (len = 0; len <= 140; len += (len < 35 || len >= 120) ? 1 : 17)
        for (off = 0; off < 8; off += (off < 2) ? 1 : 5)
        for (alias = 0; alias < 3; alias++)
        {
            nn_ptr out = (alias == 0) ? r : (alias == 1) ? a : b;

            for (i = 0; i < MUL_MAXLEN + 16; i++)
                a[i] = b[i] = r[i] = guard;
            for (i = 0; i < len; i++)
            {
                a0[i] = a[off + i] = _mul_residue(state, n);
                b0[i] = b[off + i] = _mul_residue(state, n);
            }

            _nmod_vec_mul(out + off, a + off, b + off, len, mod);
            if (!_mul_check(out, a0, b0, off, len, mod, guard))
                TEST_FUNCTION_FAIL("nu = %d, n = %wu, len = %wd, off = %wd, alias = %d\n",
                                   nu, n, len, off, alias);
        }
    }

    /* random moduli and lengths, vectors at different offsets */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        nmod_t mod;
        ulong n = n_randtest_not_zero(state);
        slong offa = n_randint(state, 8), offb = n_randint(state, 8), offr = n_randint(state, 8);

        len = n_randint(state, MUL_MAXLEN);
        nmod_init(&mod, n);

        for (i = 0; i < MUL_MAXLEN + 16; i++)
            a[i] = b[i] = r[i] = guard;
        for (i = 0; i < len; i++)
        {
            a0[i] = a[offa + i] = _mul_residue(state, n);
            b0[i] = b[offb + i] = _mul_residue(state, n);
        }

        _nmod_vec_mul(r + offr, a + offa, b + offb, len, mod);
        if (!_mul_check(r, a0, b0, offr, len, mod, guard))
            TEST_FUNCTION_FAIL("n = %wu, len = %wd, offsets %wd, %wd, %wd\n",
                               n, len, offa, offb, offr);
    }

    flint_aligned_free(a);
    flint_aligned_free(b);
    flint_aligned_free(r);
    flint_free(a0);
    flint_free(b0);

    TEST_FUNCTION_END(state);
}
