/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "nmod.h"
#include "nmod_vec.h"

TEST_FUNCTION_START(nmod_vec_dot, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong len;
        nmod_t mod;
        ulong m, res;
        nn_ptr x, y;
        mpz_t s, t;
        slong j;

        len = n_randint(state, 1000) + 1;
        m = n_randtest_not_zero(state);

        nmod_init(&mod, m);

        x = _nmod_vec_init(len);
        y = _nmod_vec_init(len);

        _nmod_vec_randtest(x, state, len, mod);
        _nmod_vec_randtest(y, state, len, mod);

        const dot_params_t params = _nmod_vec_dot_params(len, mod);

        res = _nmod_vec_dot(x, y, len, mod, params);

        mpz_init(s);
        mpz_init(t);

        for (j = 0; j < len; j++)
        {
            flint_mpz_set_ui(t, x[j]);
            flint_mpz_addmul_ui(s, t, y[j]);
        }

        flint_mpz_mod_ui(s, s, m);

        if (flint_mpz_get_ui(s) != res)
            TEST_FUNCTION_FAIL(
                    "m = %wu\n"
                    "len = %wd\n"
                    "limbs1 = %d\n",
                    m, len, params);

        mpz_clear(s);
        mpz_clear(t);

        _nmod_vec_clear(x);
        _nmod_vec_clear(y);
    }

    /*
        Two-limb sums at the largest value, all entries n - 1, at the lengths
        where the final reduction changes (nu + bitcount(len - 1) =
        FLINT_BITS - 2 and FLINT_BITS, see _NMOD_VEC_DOT2_RED) and around them;
        the expected value is len * (n-1)^2 = len mod n. Through the short
        inline code and the kernels (_nmod_vec_dot, _nmod_vec_dot_rev) for
        short lengths, and through NMOD_VEC_DOT (no memory needed) for all.
    */
    {
        flint_bitcnt_t nu;
        nn_ptr x = _nmod_vec_init(32);

        for (nu = FLINT_BITS / 2 + 1; nu <= FLINT_BITS; nu++)
        {
            int k;

            for (k = 0; k < 3; k++)
            {
                nmod_t mod;
                ulong n, res;
                slong len, j, t;
                slong lens[24];
                int nlens = 0;

                /* n of nu bits, not a power of two */
                if (k == 0)
                    n = (UWORD(1) << (nu - 1)) + 1;
                else if (k == 1)
                    n = (nu == FLINT_BITS) ? UWORD_MAX : (UWORD(1) << nu) - 1;
                else
                    n = n_randbits(state, nu) | 1;

                nmod_init(&mod, n);

                for (len = 1; len <= 12; len++)
                    lens[nlens++] = len;
                for (t = FLINT_BITS - 2; t <= FLINT_BITS; t += 2)
                {
                    /* the largest length with nu + bitcount(len - 1) <= t,
                       and the next two */
                    if ((slong) t - (slong) nu >= 4 && (slong) t - (slong) nu <= 16)
                    {
                        len = WORD(1) << (t - nu);
                        lens[nlens++] = len;
                        lens[nlens++] = len + 1;
                        lens[nlens++] = len + 2;
                    }
                }

                for (j = 0; j < 32; j++)
                    x[j] = n - 1;

                for (t = 0; t < nlens; t++)
                {
                    volatile ulong a = n - 1;
                    slong i;
                    dot_params_t params;

                    len = lens[t];
                    params = _nmod_vec_dot_params(len, mod);

                    NMOD_VEC_DOT(res, i, len, a, a, mod, params);
                    if (res != (ulong) len % n)
                        TEST_FUNCTION_FAIL("NMOD_VEC_DOT, all entries n - 1\n"
                                "n = %wu, len = %wd, method %d: got %wu, expected %wu\n",
                                n, len, (int) params.method, res, (ulong) len % n);

                    if (len <= 32)
                    {
                        res = _nmod_vec_dot(x, x, len, mod, params);
                        if (res != (ulong) len % n)
                            TEST_FUNCTION_FAIL("_nmod_vec_dot, all entries n - 1\n"
                                    "n = %wu, len = %wd, method %d: got %wu, expected %wu\n",
                                    n, len, (int) params.method, res, (ulong) len % n);

                        res = _nmod_vec_dot_rev(x, x, len, mod, params);
                        if (res != (ulong) len % n)
                            TEST_FUNCTION_FAIL("_nmod_vec_dot_rev, all entries n - 1\n"
                                    "n = %wu, len = %wd, method %d: got %wu, expected %wu\n",
                                    n, len, (int) params.method, res, (ulong) len % n);
                    }
                }
            }
        }

        _nmod_vec_clear(x);
    }

    /*
        The two-limb code called directly for n <= 3, which
        _nmod_vec_dot_params never does (one limb suffices there): the
        kernels _nmod_vec_dot2 and _nmod_vec_dot2_rev, and the short inline
        code through _nmod_vec_dot with the method _DOT2
    */
    {
        ulong n;
        nn_ptr x = _nmod_vec_init(20);
        nn_ptr y = _nmod_vec_init(20);

        for (n = 1; n <= 3; n++)
        {
            nmod_t mod;
            slong len, j;

            nmod_init(&mod, n);

            for (len = 1; len <= 20; len++)
            {
                ulong res, s1 = 0, s2 = 0;
                dot_params_t params = {_DOT2, UWORD(0)};

                for (j = 0; j < len; j++)
                {
                    x[j] = n_randint(state, n);
                    y[j] = n_randint(state, n);
                }
                for (j = 0; j < len; j++)
                {
                    s1 += x[j] * y[j];
                    s2 += x[j] * y[len - 1 - j];
                }
                s1 %= n;
                s2 %= n;

                res = _nmod_vec_dot2(x, y, len, mod);
                if (res != s1)
                    TEST_FUNCTION_FAIL("_nmod_vec_dot2, n = %wu, len = %wd: got %wu, expected %wu\n",
                                       n, len, res, s1);

                res = _nmod_vec_dot2_rev(x, y, len, mod);
                if (res != s2)
                    TEST_FUNCTION_FAIL("_nmod_vec_dot2_rev, n = %wu, len = %wd: got %wu, expected %wu\n",
                                       n, len, res, s2);

                if (len <= 11)
                {
                    res = _nmod_vec_dot(x, y, len, mod, params);
                    if (res != s1)
                        TEST_FUNCTION_FAIL("_nmod_vec_dot with _DOT2, n = %wu, len = %wd: got %wu, expected %wu\n",
                                           n, len, res, s1);
                }
            }
        }

        _nmod_vec_clear(x);
        _nmod_vec_clear(y);
    }

#if FLINT_BITS == 64
    /*
        One long dot product through the split code (_DOT2_SPLIT), which has
        no length bound: all entries n - 1, with n close to the largest
        modulus of _DOT2_SPLIT and 2^56 mod n large, through the NMOD_VEC_DOT
        macro (no memory needed); the expected value is
        len * (n-1)^2 = len mod n. Computing the final value
        pow2_precomp * hi + lo in a single word would overflow here. The
        split code is used even where _nmod_vec_dot_params does not choose it
        for this length (AArch64 uses _DOT2 from length
        NMOD_VEC_DOT2_SPLIT_MAX_LEN on).
    */
    {
        nmod_t mod;
        dot_params_t params;
        volatile ulong a;
        ulong res, n = UWORD(1515531501);
        slong k, len = WORD(300000007);

        nmod_init(&mod, n);
#if !defined(NMOD_VEC_DOT2_SPLIT_MAX_LEN)
        params = _nmod_vec_dot_params(len, mod);
        if (params.method != _DOT2_SPLIT)
            TEST_FUNCTION_FAIL("n = %wu, len = %wd: method %d, expected _DOT2_SPLIT\n",
                               n, len, (int) params.method);
#endif
        params.method = _DOT2_SPLIT;
        NMOD_RED(params.pow2_precomp, UWORD(1) << DOT_SPLIT_BITS, mod);

        a = n - 1;
        NMOD_VEC_DOT(res, k, len, a, a, mod, params);

        if (res != (ulong) len % n)
            TEST_FUNCTION_FAIL("long split dot product, n = %wu, len = %wd: got %wu, expected %wu\n",
                               n, len, res, (ulong) len % n);
    }
#endif

    TEST_FUNCTION_END(state);
}
