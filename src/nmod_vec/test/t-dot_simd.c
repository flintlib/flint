/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "ulong_extras.h"
#include "nmod_vec.h"

/*
    The SIMD dot products for moduli above 2^32 (_nmod_vec_dot_u52,
    _nmod_vec_dot_split_limbs and _nmod_vec_dot_u64, forward / reversed / through
    row pointers), called directly on their whole range of moduli, entries
    at the maximum, and lengths beyond the internal chunks of the
    accumulators; on machines where they are not available they fall back
    to the scalar code and the test still runs. Also checks that
    _nmod_vec_dot_params returns them where expected.
*/

static ulong
dot_ref(nn_srcptr x, nn_srcptr y, slong len, ulong m, int rev)
{
    mpz_t s, t;
    slong j;
    ulong r;

    mpz_init(s);
    mpz_init(t);
    for (j = 0; j < len; j++)
    {
        flint_mpz_set_ui(t, x[j]);
        flint_mpz_addmul_ui(s, t, rev ? y[len - 1 - j] : y[j]);
    }
    flint_mpz_mod_ui(s, s, m);
    r = flint_mpz_get_ui(s);
    mpz_clear(s);
    mpz_clear(t);

    return r;
}

TEST_FUNCTION_START(nmod_vec_dot_simd, state)
{
    slong i;

    for (i = 0; i < 300 * flint_test_multiplier(); i++)
    {
        slong len, j, offset;
        ulong m, r0, r1, a, b, c;
        nmod_t mod;
        nn_ptr x, y, store;
        nn_ptr * rows;
        dot_params_t params;
        int u52, split;

        switch (n_randint(state, 9))
        {
            case 0:
                m = n_randbits(state, 33 + n_randint(state, 20));
                break;
            case 1:
                m = (UWORD(1) << (32 + n_randint(state, 20))) + 1 + n_randint(state, 4);
                break;
            case 2:
                m = (UWORD(1) << (33 + n_randint(state, 20))) - 1 - n_randint(state, 4);
                break;
            case 3:
                m = UWORD(1) << 52;
                break;
            case 4:
                m = n_randbits(state, 53 + n_randint(state, 12));
                break;
            case 5:
                m = (UWORD(1) << 52) + 1 + n_randint(state, 4);
                break;
            case 6:
                /* around the limits of _DOT_SPLIT_LIMBS with IFMA (58 bits)
                   and of its two variants (62 bits) */
                m = (UWORD(1) << (58 + 4 * n_randint(state, 2)))
                        + n_randint(state, 5) - 2;
                break;
            case 7:
                m = (UWORD(1) << (60 + n_randint(state, 4))) + n_randint(state, 5) - 2;
                break;
            default:
                m = UWORD_MAX - n_randint(state, 4);
                break;
        }
        if (m < 2)
            m = 2;
        u52 = (m <= (UWORD(1) << 52));
        split = NMOD_VEC_HAVE_DOT_SPLIT_LIMBS
                && (!NMOD_VEC_HAVE_DOT_U64 || m <= (UWORD(1) << 58));

        switch (n_randint(state, 4))
        {
            case 0: len = n_randint(state, 64); break;
            case 1: len = n_randint(state, 1000); break;
            case 2: len = n_randint(state, 20000); break;
            default: len = 131072 - 40 + n_randint(state, 80); break;
        }

        nmod_init(&mod, m);
        x = _nmod_vec_init(len + 1);
        y = _nmod_vec_init(len + 1);
        offset = n_randint(state, 3);
        store = _nmod_vec_init(3 * (len + 1));
        rows = flint_malloc((len + 1) * sizeof(nn_ptr));

        switch (n_randint(state, 3))
        {
            case 0:
                _nmod_vec_randtest(x, state, len, mod);
                _nmod_vec_randtest(y, state, len, mod);
                break;
            case 1:
                for (j = 0; j < len; j++)
                    x[j] = y[j] = m - 1;
                break;
            default:
                for (j = 0; j < len; j++)
                {
                    x[j] = m - 1 - n_randint(state, 2);
                    y[j] = m - 1 - n_randint(state, 2);
                }
                break;
        }
        for (j = 0; j < len; j++)
        {
            rows[j] = store + 3 * j;
            rows[j][offset] = y[j];
        }

        r0 = dot_ref(x, y, len, m, 0);
        r1 = dot_ref(x, y, len, m, 1);

        if (u52)
        {
            a = _nmod_vec_dot_u52(x, y, len, mod);
            b = _nmod_vec_dot_u52_rev(x, y, len, mod);
            c = _nmod_vec_dot_u52_ptr(x, rows, offset, len, mod);
            if (a != r0 || b != r1 || c != r0)
                TEST_FUNCTION_FAIL("u52: m = %wu, len = %wd: %wu %wu %wu, "
                                   "expected %wu %wu\n", m, len, a, b, c, r0, r1);
        }

        a = _nmod_vec_dot_u64(x, y, len, mod);
        b = _nmod_vec_dot_u64_rev(x, y, len, mod);
        c = _nmod_vec_dot_u64_ptr(x, rows, offset, len, mod);
        if (a != r0 || b != r1 || c != r0)
            TEST_FUNCTION_FAIL("u64: m = %wu, len = %wd: %wu %wu %wu, "
                               "expected %wu %wu\n", m, len, a, b, c, r0, r1);

        a = _nmod_vec_dot_split_limbs(x, y, len, mod);
        b = _nmod_vec_dot_split_limbs_rev(x, y, len, mod);
        c = _nmod_vec_dot_split_limbs_ptr(x, rows, offset, len, mod);
        if (a != r0 || b != r1 || c != r0)
            TEST_FUNCTION_FAIL("split_limbs: m = %wu, len = %wd: %wu %wu %wu, "
                               "expected %wu %wu\n", m, len, a, b, c, r0, r1);

        /* the dispatch, whatever it chooses */
        params = _nmod_vec_dot_params(len, mod);
        a = _nmod_vec_dot(x, y, len, mod, params);
        b = _nmod_vec_dot_rev(x, y, len, mod, params);
        c = _nmod_vec_dot_ptr(x, rows, offset, len, mod, params);
        if (a != r0 || b != r1 || c != r0)
            TEST_FUNCTION_FAIL("dispatch (method %d): m = %wu, len = %wd\n",
                               params.method, m, len);

        if ((m & (m - 1)) != 0 && _nmod_vec_dot_bound_limbs(len, mod) <= 2)
        {
            dot_method_t exp = _DOT2;
            if (NMOD_VEC_HAVE_DOT_U52 && u52 && len >= NMOD_VEC_DOT_U52_MIN_LEN)
                exp = _DOT_U52;
            else if (split && len >= NMOD_VEC_DOT_SPLIT_LIMBS_MIN_LEN)
                exp = _DOT_SPLIT_LIMBS;
            else if (NMOD_VEC_HAVE_DOT_U64 && !u52 && len >= NMOD_VEC_DOT_U64_MIN_LEN)
                exp = _DOT_U64;
            if (exp != _DOT2 && params.method != exp)
                TEST_FUNCTION_FAIL("params: expected %d, got %d: "
                                   "m = %wu, len = %wd\n", exp, params.method, m, len);
        }
        else if ((m & (m - 1)) != 0)
        {
            dot_method_t exp = _DOT3;
            if (split && len >= NMOD_VEC_DOT_SPLIT_LIMBS_MIN_LEN)
                exp = _DOT3_SPLIT_LIMBS;
            else if (NMOD_VEC_HAVE_DOT_U64 && len >= NMOD_VEC_DOT_U64_MIN_LEN)
                exp = _DOT3_U64;
            if (exp != _DOT3 && params.method != exp)
                TEST_FUNCTION_FAIL("params: expected %d, got %d: "
                                   "m = %wu, len = %wd\n", exp, params.method, m, len);
        }

        _nmod_vec_clear(x);
        _nmod_vec_clear(y);
        _nmod_vec_clear(store);
        flint_free(rows);
    }

    TEST_FUNCTION_END(state);
}
