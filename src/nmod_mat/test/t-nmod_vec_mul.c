/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_vec.h"
#include "nmod_mat.h"

TEST_FUNCTION_START(nmod_mat_nmod_vec_mul, state)
{
    slong i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        ulong p;
        nmod_mat_t A, B, C;
        ulong * a, * c;
        ulong ** aa, ** cc;
        slong j, m, n, alen;

        p = n_randtest_not_zero(state);
        m = n_randint(state, 50);
        n = n_randint(state, 50);
        alen = n_randint(state, 50);

        nmod_mat_init(C, 1, n, p);
        nmod_mat_init(A, 1, m, p);
        nmod_mat_init(B, m, n, p);
        c = _nmod_vec_init(n);
        a = _nmod_vec_init(alen);

        nmod_mat_randtest(B, state);
        _nmod_vec_randtest(c, state, n, B->mod);
        _nmod_vec_randtest(a, state, alen, B->mod);

        cc = FLINT_ARRAY_ALLOC(n, ulong*);
        for (j = 0; j < n; j++)
        {
            cc[j] = FLINT_ARRAY_ALLOC(1, ulong);
            cc[j][0] = c[j];
        }

        aa = FLINT_ARRAY_ALLOC(alen, ulong*);
        for (j = 0; j < alen; j++)
        {
            aa[j] = FLINT_ARRAY_ALLOC(1, ulong);
            aa[j][0] = a[j];
        }

        nmod_mat_nmod_vec_mul(c, a, alen, B);
        nmod_mat_nmod_vec_mul_ptr(cc, (const ulong * const *)aa, alen, B);

        /* supposed to match mul of the chopped or zero-extended a */
        for (j = 0; j < m && j < alen; j++)
            nmod_mat_entry(A, 0, j) = a[j];

        nmod_mat_mul(C, A, B);

        for (j = 0; j < n; j++)
            if (nmod_mat_entry(C, 0, j) != c[j] || nmod_mat_entry(C, 0, j) != cc[j][0])
                TEST_FUNCTION_FAIL("");

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        _nmod_vec_clear(c);
        _nmod_vec_clear(a);

        for (j = 0; j < n; j++)
        {
            flint_free(cc[j]);
        }
        flint_free(cc);

        for (j = 0; j < alen; j++)
        {
            flint_free(aa[j]);
        }
        flint_free(aa);
    }

    /* against one modular multiplication per entry, on the dimensions
       and moduli of the vectorized paths: rows beyond the chunks of the
       accumulators, columns beyond the blocks, moduli at the limits of
       each method, entries n - 1 (largest products) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_mat_t B;
        ulong * a, * c, * d;
        ulong n;
        slong j, k, m;
        int mode;

        switch (n_randint(state, 8))
        {
            case 0: n = UWORD(1) << 26; break;
            case 1: n = (UWORD(1) << 26) + 1; break;
            case 2: n = UWORD(1) << 32; break;
            case 3: n = (UWORD(1) << 50) - 1; break;
            case 4: n = UWORD(1) << 52; break;
            case 5: n = (UWORD(1) << 52) + 1; break;
            default: n = n_randbits(state, 2 + n_randint(state, FLINT_BITS - 1)); break;
        }
        n = FLINT_MAX(n, UWORD(2));

        k = 1 + n_randint(state, n_randint(state, 4) ? 80 : 1200);
        m = 1 + n_randint(state, n_randint(state, 4) ? 80 : 4500);
        mode = n_randint(state, 3);

        nmod_mat_init(B, k, m, n);
        a = _nmod_vec_init(k);
        c = _nmod_vec_init(m + 1);
        d = _nmod_vec_init(m);

        for (j = 0; j < k * m; j++)
            B->entries[j] = mode ? n - 1 - (mode == 2 ? n_randint(state, 2) : 0)
                                 : n_randint(state, n);
        for (j = 0; j < k; j++)
            a[j] = mode ? n - 1 - (mode == 2 ? n_randint(state, 2) : 0)
                        : n_randint(state, n);
        c[m] = 12345;

        nmod_mat_nmod_vec_mul(c, a, k, B);

        _nmod_vec_scalar_mul_nmod(d, nmod_mat_entry_ptr(B, 0, 0), m, a[0], B->mod);
        for (j = 1; j < k; j++)
            _nmod_vec_scalar_addmul_nmod(d, nmod_mat_entry_ptr(B, j, 0), m, a[j], B->mod);

        if (!_nmod_vec_equal(c, d, m) || c[m] != 12345)
            TEST_FUNCTION_FAIL("k: %wd, m: %wd, n: %wu, mode: %d\n", k, m, n, mode);

        nmod_mat_clear(B);
        _nmod_vec_clear(a);
        _nmod_vec_clear(c);
        _nmod_vec_clear(d);
    }

    TEST_FUNCTION_END(state);
}
