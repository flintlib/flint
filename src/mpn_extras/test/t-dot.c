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
#include "fmpz.h"

/* {res, s} = sum a[k] b[k] mod 2^(64 s) using fmpz arithmetic */
static void
_dot_ref(nn_ptr res, nn_srcptr a, slong n1, nn_srcptr b, slong n2, slong len, slong s, int sgn)
{
    fmpz_t x, y, acc;
    slong k;

    fmpz_init(x);
    fmpz_init(y);
    fmpz_init(acc);

    for (k = 0; k < len; k++)
    {
        if (sgn)
        {
            fmpz_set_signed_ui_array(x, a + k * n1, n1);
            fmpz_set_signed_ui_array(y, b + k * n2, n2);
        }
        else
        {
            fmpz_set_ui_array(x, a + k * n1, n1);
            fmpz_set_ui_array(y, b + k * n2, n2);
        }

        fmpz_addmul(acc, x, y);
    }

    fmpz_fdiv_r_2exp(acc, acc, s * FLINT_BITS);
    fmpz_get_ui_array(res, s, acc);

    fmpz_clear(x);
    fmpz_clear(y);
    fmpz_clear(acc);
}

TEST_FUNCTION_START(flint_mpn_dot, state)
{
    slong iter;

    for (iter = 0; iter < 20000 * flint_test_multiplier(); iter++)
    {
        slong n1, n2, s, len, i, astride, bstride, k;
        nn_ptr a, b, ap, bp, r1, r2, t;
        int sgn, kind;

        n1 = 1 + n_randint(state, 6);
        n2 = 1 + n_randint(state, n1);
        s = n1 + n2 - 1 + n_randint(state, 3);
        len = 1 + n_randint(state, 12);
        if (n_randint(state, 20) == 0)
            len = 1 + n_randint(state, 200);
        sgn = n_randint(state, 2);
        /* 0 = forward, 1 = reversed, 2 = strided, 3 = generic */
        kind = n_randint(state, 4);

        /* strides in units of elements */
        astride = (kind == 2) ? 1 + n_randint(state, 3) : 1;
        bstride = (kind == 2) ? 1 + n_randint(state, 3) : 1;

        a = flint_malloc(sizeof(ulong) * n1 * len * astride);
        b = flint_malloc(sizeof(ulong) * n2 * len * bstride);
        ap = flint_malloc(sizeof(ulong) * n1 * len);
        bp = flint_malloc(sizeof(ulong) * n2 * len);
        r1 = flint_malloc(sizeof(ulong) * s);
        r2 = flint_malloc(sizeof(ulong) * s);
        t = flint_malloc(sizeof(ulong) * (n1 + n2));

        for (i = 0; i < n1 * len * astride; i++)
            a[i] = n_randtest(state);
        for (i = 0; i < n2 * len * bstride; i++)
            b[i] = n_randtest(state);
        if (n_randint(state, 2))
            flint_mpn_rrandom(a, state, n1 * len * astride);
        if (n_randint(state, 2))
            flint_mpn_rrandom(b, state, n2 * len * bstride);
        for (i = 0; i < s; i++)
            r1[i] = n_randtest(state);

        /* packed copies in the order of the products, for the reference */
        for (k = 0; k < len; k++)
        {
            flint_mpn_copyi(ap + k * n1, a + k * astride * n1, n1);
            if (kind == 1)
                flint_mpn_copyi(bp + k * n2, b + (len - 1 - k) * n2, n2);
            else
                flint_mpn_copyi(bp + k * n2, b + k * bstride * n2, n2);
        }

        {
            flint_mpn_dot_func_t fdot = NULL;
            flint_mpn_dot_strided_func_t fstr = NULL;

            if (s - (n1 + n2 - 1) < 3)
            {
                if (kind == 0 && n1 <= FLINT_MPN_DOT_DEDICATED_TAB_N)
                    fdot = flint_mpn_dot_tab[sgn][n1][n2][s - (n1 + n2 - 1)];
                else if (kind == 1 && n1 <= FLINT_MPN_DOT_DEDICATED_TAB_N)
                    fdot = flint_mpn_dot_rev_tab[sgn][n1][n2][s - (n1 + n2 - 1)];
                else if (kind <= 2 && n1 <= FLINT_MPN_DOT_TAB_N)
                    fstr = flint_mpn_dot_strided_tab[sgn][n1][n2][s - (n1 + n2 - 1)];
            }

            if (fdot)
                fdot(r1, a, b, len);
            else if (fstr && kind == 1)
            {
                /* reversed access through the strided kernels, as in the
                   dot_rev dispatchers for the larger fixed sizes */
                fstr(r1, a, n1, b + (len - 1) * n2, -n2, len);
            }
            else if (fstr)
                fstr(r1, a, astride * n1, b, bstride * n2, len);
        else
        {
            /* generic (reversed): use the packed copies, reversing b */
            nn_ptr brev = flint_malloc(sizeof(ulong) * n2 * len);
            for (k = 0; k < len; k++)
                flint_mpn_copyi(brev + k * n2, bp + (len - 1 - k) * n2, n2);
            if (sgn)
                _flint_mpn_dot_rev_generic_signed(r1, ap, n1, brev, n2, len, s, t);
            else
                _flint_mpn_dot_rev_generic(r1, ap, n1, brev, n2, len, s, t);
            flint_free(brev);
        }
        }

        _dot_ref(r2, ap, n1, bp, n2, len, s, sgn);

        if (mpn_cmp(r1, r2, s) != 0)
            TEST_FUNCTION_FAIL(
                    "n1 = %wd, n2 = %wd, s = %wd, len = %wd, sgn = %d, kind = %d, astride = %wd, bstride = %wd\n"
                    "ap = %{ulong*}\n"
                    "bp = %{ulong*}\n"
                    "r1 = %{ulong*}\n"
                    "r2 = %{ulong*}\n",
                    n1, n2, s, len, sgn, kind, astride, bstride, ap, n1 * len, bp, n2 * len, r1, s, r2, s);

        flint_free(a);
        flint_free(b);
        flint_free(ap);
        flint_free(bp);
        flint_free(r1);
        flint_free(r2);
        flint_free(t);
    }

    TEST_FUNCTION_END(state);
}
