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

TEST_FUNCTION_START(nfixed_mat_mul, state)
{
    slong iter, m, n, p, i, nlimbs;
    nn_ptr A, B, C, D, t;
    nn_ptr a;
    int which;

    slong MAXN = 12;
    slong MINLIMBS = 2;
    slong MAXLIMBS = 12;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        which = n_randint(state, 6);

        m = 1 + n_randint(state, MAXN);
        n = 1 + n_randint(state, MAXN);
        p = 1 + n_randint(state, MAXN);

        nlimbs = MINLIMBS + n_randint(state, MAXLIMBS - MINLIMBS + 1);

        ulong maxerr = 2 * (2 * nlimbs - 1) * n;

        A = flint_malloc((nlimbs + 1) * (m * n) * sizeof(ulong));
        B = flint_malloc((nlimbs + 1) * (n * p) * sizeof(ulong));
        C = flint_malloc((nlimbs + 1) * (m * p) * sizeof(ulong));
        D = flint_malloc((nlimbs + 1) * (m * p) * sizeof(ulong));
        t = flint_malloc((nlimbs + 1) * sizeof(ulong));

        for (i = 0; i < m * n; i++)
        {
            a = A + i * (nlimbs + 1);
            a[0] = n_randint(state, 2);
            flint_mpn_rrandom(a + 1, state, nlimbs);
            a[nlimbs] >>= 10;
        }

        for (i = 0; i < n * p; i++)
        {
            a = B + i * (nlimbs + 1);
            a[0] = n_randint(state, 2);
            flint_mpn_rrandom(a + 1, state, nlimbs);
            a[nlimbs] >>= 10;
        }

        for (i = 0; i < m * p; i++)
        {
            a = C + i * (nlimbs + 1);
            a[0] = n_randint(state, 2);
            flint_mpn_rrandom(a + 1, state, nlimbs);

            a = D + i * (nlimbs + 1);
            a[0] = n_randint(state, 2);
            flint_mpn_rrandom(a + 1, state, nlimbs);
        }

        _nfixed_mat_mul_classical(C, A, B, m, n, p, nlimbs);

        if (which == 0)
            _nfixed_mat_mul_waksman(D, A, B, m, n, p, nlimbs);
        else
            _nfixed_mat_mul_strassen(D, A, B, m, n, p, which, nlimbs);

        for (i = 0; i < m * p; i++)
        {
            nfixed_sub(t, C + i * (nlimbs + 1), D + i * (nlimbs + 1), nlimbs);

            if (!flint_mpn_zero_p(t + 2, nlimbs - 1) || t[1] > maxerr)
            {
                TEST_FUNCTION_FAIL("nlimbs = %wd, m = %wd, n = %wd, p = %wd\n\nt = %{ulong*}, maxerr = %wu\n\nA = %{ulong*}\n\nB = %{ulong*}\n\nC = %{ulong*}\n\nD = %{ulong*}\n\n",
                    nlimbs, m, n, p,
                    t, nlimbs + 1, maxerr, A, m * n * (nlimbs + 1), B, n * p * (nlimbs + 1), C, m * p * (nlimbs + 1), D, m * p * (nlimbs + 1));
            }
        }

        flint_free(A);
        flint_free(B);
        flint_free(C);
        flint_free(D);
    }

    TEST_FUNCTION_END(state);
}