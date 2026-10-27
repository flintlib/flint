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
