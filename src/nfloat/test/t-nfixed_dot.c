/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "arf.h"
#include "gr_vec.h"
#include "gr_special.h"
#include "nfloat.h"

#define MAXLEN 10
#define MINLIMBS 2
#define MAXLIMBS 8

FLINT_FORCE_INLINE
void nfixed_add(nn_ptr res, nn_srcptr a, nn_srcptr b, slong nlimbs)
{
    int asgn, bsgn;
    asgn = a[0];
    bsgn = b[0];

    if (asgn == bsgn)
    {
        res[0] = asgn;
        mpn_add_n(res + 1, a + 1, b + 1, nlimbs);
    }
    else
    {
        res[0] = asgn ^ flint_mpn_signed_sub_n(res + 1, a + 1, b + 1, nlimbs);
    }
}

FLINT_FORCE_INLINE
void nfixed_sub(nn_ptr res, nn_srcptr a, nn_srcptr b, slong nlimbs)
{
    int asgn, bsgn;
    asgn = a[0];
    bsgn = b[0];

    if (asgn != bsgn)
    {
        res[0] = asgn;
        mpn_add_n(res + 1, a + 1, b + 1, nlimbs);
    }
    else
    {
        res[0] = asgn ^ flint_mpn_signed_sub_n(res + 1, a + 1, b + 1, nlimbs);
    }
}

FLINT_FORCE_INLINE
void nfixed_mul(nn_ptr res, nn_srcptr a, nn_srcptr b, slong nlimbs)
{
    int asgn, bsgn;
    asgn = a[0];
    bsgn = b[0];
    res[0] = asgn ^ bsgn;
    flint_mpn_mulhigh_n(res + 1, a + 1, b + 1, nlimbs);
}

TEST_FUNCTION_START(nfixed_dot, state)
{
    slong iter, len, i, nlimbs;
    nn_ptr a;

    ulong A[MAXLEN * (MAXLIMBS + 1)];
    ulong B[MAXLEN * (MAXLIMBS + 1)];
    ulong C[MAXLIMBS + 1];
    ulong D[MAXLIMBS + 1];
    ulong t[MAXLIMBS + 1];

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        len = 1 + n_randint(state, MAXLEN);
        nlimbs = MINLIMBS + n_randint(state, MAXLIMBS - MINLIMBS + 1);

        ulong maxerr = (2 * nlimbs - 1) * len;

        for (i = 0; i < len; i++)
        {
            a = A + i * (nlimbs + 1);
            a[0] = n_randint(state, 2);
            flint_mpn_rrandom(a + 1, state, nlimbs);
            a[nlimbs] >>= 10;

            a = B + i * (nlimbs + 1);
            a[0] = n_randint(state, 2);
            flint_mpn_rrandom(a + 1, state, nlimbs);
            a[nlimbs] >>= 10;
        }

        switch (nlimbs)
        {
            case 2:
                _nfixed_dot_2(C, A, nlimbs + 1, B, nlimbs + 1, len);
                break;
            case 3:
                _nfixed_dot_3(C, A, nlimbs + 1, B, nlimbs + 1, len);
                break;
            case 4:
                _nfixed_dot_4(C, A, nlimbs + 1, B, nlimbs + 1, len);
                break;
            case 5:
                _nfixed_dot_5(C, A, nlimbs + 1, B, nlimbs + 1, len);
                break;
            case 6:
                _nfixed_dot_6(C, A, nlimbs + 1, B, nlimbs + 1, len);
                break;
            case 7:
                _nfixed_dot_7(C, A, nlimbs + 1, B, nlimbs + 1, len);
                break;
            case 8:
                _nfixed_dot_8(C, A, nlimbs + 1, B, nlimbs + 1, len);
                break;
            default:
                flint_abort();
        }

        flint_mpn_zero(D, nlimbs + 1);

        for (i = 0; i < len; i++)
        {
            nfixed_mul(t, A + i * (nlimbs + 1), B + i * (nlimbs + 1), nlimbs);
            nfixed_add(D, D, t, nlimbs);
        }

        nfixed_sub(t, C, D, nlimbs);
        if (!flint_mpn_zero_p(t + 2, nlimbs - 1) || t[1] > maxerr)
        {
            TEST_FUNCTION_FAIL("nlimbs = %wd, len = %wd,\n\nA = %{ulong*},\n\nB = %{ulong*},\n\nC = %{ulong*},\n\nD = %{ulong*},\n\nt = %{ulong*}\n", nlimbs, len,
                A, len * (nlimbs + 1), B, len * (nlimbs + 1),
                C, nlimbs + 1,
                D, nlimbs + 1,
                t, nlimbs + 1);
        }
    }

    TEST_FUNCTION_END(state);
}

#undef MAXLEN
#undef MINLIMBS
#undef MAXLIMBS
