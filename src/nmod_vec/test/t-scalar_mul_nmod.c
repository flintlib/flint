/*
    Copyright (C) 2009 William Hart

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

/* a residue, often at the corners: n - 1 - d for small d, small, or random */
static ulong
_smul_residue(flint_rand_t state, ulong n)
{
    ulong d;

    switch (n_randint(state, 4))
    {
        case 0:
            d = n_randint(state, 3);
            return (d < n) ? n - 1 - d : n - 1;
        case 1:
            return n_randint(state, 3) % n;
        default:
            return n_randint(state, n);
    }
}

/* the moduli 2^(nu - 1), 2^(nu - 1) + 1, 2^nu - 1 and a random one */
static ulong
_smul_modulus(flint_rand_t state, int nu, int k)
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

#define SMUL_MAXLEN 160

/* scalar_mul and scalar_addmul, and their lazy variants up to 62 bits,
   against nmod_mul and nmod_add: every bit size, c at the corners, lengths
   0 to 3 vectors of 8 words plus a tail and around the peeling threshold of
   128, offsets 0, 1, 2 and 7 words from a 64-byte boundary, in place and
   not; the lazy outputs must be below 2n, and the lazy scalar_addmul gets
   res in [0, 2n) (reduced in place, where it is also the input vec); the
   words around the output must not change; returns 0 on success */
static int
_smul_check_all(flint_rand_t state)
{
    const ulong guard = UWORD_MAX / 3;
    nn_ptr a, r, r0, want;
    slong i, len, off;
    int nu, k, ci, add, inplace, lazy, bad = 0;

    a = flint_aligned_alloc(64, (SMUL_MAXLEN + 16) * sizeof(ulong));
    r = flint_aligned_alloc(64, (SMUL_MAXLEN + 16) * sizeof(ulong));
    r0 = flint_malloc((SMUL_MAXLEN + 16) * sizeof(ulong));
    want = flint_malloc(SMUL_MAXLEN * sizeof(ulong));

    for (nu = 2; nu <= FLINT_BITS; nu++)
    for (k = 0; k < 4; k++)
    {
        nmod_t mod;
        ulong n = _smul_modulus(state, nu, k);
        nmod_init(&mod, n);

        for (ci = 0; ci < 6; ci++)
        {
            ulong c = (ci == 0) ? 2 % n : (ci == 1) ? n - 2 : (ci == 2) ? n / 2
                    : (ci == 3) ? n_randint(state, n) : (ci == 4) ? 1 : n - 1;

            for (len = 0; len <= 140; len += (len < 35 || len >= 120) ? 1 : 17)
            for (off = 0; off < 8; off += (off < 2) ? 1 : 5)
            for (add = 0; add < 2; add++)
            for (inplace = 0; inplace < 2; inplace++)
            for (lazy = 0; lazy < 1 + (nu <= 62); lazy++)
            {
                nn_srcptr src = inplace ? r + off : a + off;

                for (i = 0; i < SMUL_MAXLEN + 16; i++)
                {
                    a[i] = guard;
                    r[i] = (add && i >= off && i < off + len) ? _smul_residue(state, n) : guard;
                    if (add && lazy && !inplace && i >= off && i < off + len
                            && n_randint(state, 2))
                        r[i] += n;
                }
                for (i = 0; i < len; i++)
                {
                    a[off + i] = _smul_residue(state, n);
                    if (inplace && !add)
                        r[off + i] = a[off + i];
                }
                for (i = 0; i < len; i++)
                {
                    ulong x = inplace ? r[off + i] : a[off + i];
                    want[i] = nmod_mul(x, c, mod);
                    if (add)
                        want[i] = nmod_add(want[i], r[off + i] % n, mod);
                }
                for (i = 0; i < SMUL_MAXLEN + 16; i++)
                    r0[i] = r[i];

                if (add && lazy)
                    _nmod_vec_scalar_addmul_nmod_lazy(r + off, src, len, c, mod);
                else if (add)
                    _nmod_vec_scalar_addmul_nmod(r + off, src, len, c, mod);
                else if (lazy)
                    _nmod_vec_scalar_mul_nmod_lazy(r + off, src, len, c, mod);
                else
                    _nmod_vec_scalar_mul_nmod(r + off, src, len, c, mod);

                for (i = 0; i < SMUL_MAXLEN + 16; i++)
                {
                    int in = (i >= off && i < off + len);
                    ulong w = in ? want[i - off] : r0[i];
                    int good = (in && lazy) ? (r[i] / 2 < n && r[i] % n == w) : (r[i] == w);
                    if (!good && !bad)
                    {
                        flint_printf("%s%s: nu = %d, n = %wu, c = %wu, len = %wd, "
                                     "off = %wd, in place %d, i = %wd\n",
                                     add ? "addmul" : "mul", lazy ? "_lazy" : "", nu, n, c,
                                     len, off, inplace, i - off);
                        bad = 1;
                    }
                }
            }
        }
    }

    flint_aligned_free(a);
    flint_aligned_free(r);
    flint_free(r0);
    flint_free(want);
    return bad;
}

TEST_FUNCTION_START(nmod_vec_scalar_mul_nmod, state)
{
    int i, result;

    if (_smul_check_all(state))
        TEST_FUNCTION_FAIL("scalar_mul / scalar_addmul (lazy or not): wrong result (see above)\n");

    {
        nmod_t mod;
        nmod_init(&mod, 2);
        _nmod_vec_scalar_mul_nmod(NULL, NULL, 0, 1, mod);
    }

    /* Check (a + b)*c == a*c + b*c */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong len = n_randint(state, 100) + 1;
        ulong n = n_randtest_not_zero(state);
        ulong c = n_randint(state, n);
        nmod_t mod;

        nn_ptr vec = _nmod_vec_init(len);
        nn_ptr vec2 = _nmod_vec_init(len);
        nn_ptr vec3 = _nmod_vec_init(len);

        nmod_init(&mod, n);

        _nmod_vec_randtest(vec, state, len, mod);
        _nmod_vec_randtest(vec2, state, len, mod);

        _nmod_vec_add(vec3, vec, vec2, len, mod);
        _nmod_vec_scalar_mul_nmod(vec3, vec3, len, c, mod);

        _nmod_vec_scalar_mul_nmod(vec, vec, len, c, mod);
        _nmod_vec_scalar_mul_nmod(vec2, vec2, len, c, mod);
        _nmod_vec_add(vec, vec, vec2, len, mod);

        result = _nmod_vec_equal(vec, vec3, len);
        if (!result)
            TEST_FUNCTION_FAIL("len = %wd, n = %wd\n", len, n);

        _nmod_vec_clear(vec);
        _nmod_vec_clear(vec2);
        _nmod_vec_clear(vec3);
    }

    TEST_FUNCTION_END(state);
}
