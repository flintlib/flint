/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mat.h"

/* whether B is the transpose of A, entrywise */
static int
is_transpose(const nmod_mat_t B, const nmod_mat_t A)
{
    slong i, j;

    if (B->r != A->c || B->c != A->r)
        return 0;

    for (i = 0; i < B->r; i++)
        for (j = 0; j < B->c; j++)
            if (nmod_mat_entry(B, i, j) != nmod_mat_entry(A, j, i))
                return 0;

    return 1;
}

/* dimension: mostly small, sometimes beyond the blocks of transpose.c */
static slong
rand_dim(flint_rand_t state)
{
    return n_randint(state, n_randint(state, 4) ? 20 : 150);
}

TEST_FUNCTION_START(nmod_mat_transpose, state)
{
    slong m, n, mod, mod2, rep;

    /* Rectangular transpose, same modulus */
    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A, B, C;

        m = rand_dim(state);
        n = rand_dim(state);

        mod = n_randtest_not_zero(state);

        nmod_mat_init(A, m, n, mod);
        nmod_mat_init(B, n, m, mod);
        nmod_mat_init(C, m, n, mod);

        nmod_mat_randtest(A, state);
        nmod_mat_randtest(B, state);

        nmod_mat_transpose(B, A);
        nmod_mat_transpose(C, B);

        if (!is_transpose(B, A))
            TEST_FUNCTION_FAIL("B != A^T\n");

        if (!nmod_mat_equal(C, A))
            TEST_FUNCTION_FAIL("C != A\n");

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
    }

    /* Rectangular transpose, different modulus */
    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A, AT, B, BT, AT2;

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        mod = n_randtest_not_zero(state);
        mod2 = n_randtest_not_zero(state);

        nmod_mat_init(A, m, n, mod);
        nmod_mat_init(AT, n, m, mod);
        nmod_mat_init(B, m, n, mod2);
        nmod_mat_init(BT, n, m, mod2);
        nmod_mat_init(AT2, n, m, mod2);

        nmod_mat_randtest(A, state);
        nmod_mat_set(B, A);

        nmod_mat_transpose(AT, A);
        nmod_mat_transpose(BT, B);

        nmod_mat_set(AT2, AT);

        if (!nmod_mat_equal(BT, AT2))
            TEST_FUNCTION_FAIL("AT != BT\n");

        nmod_mat_clear(A);
        nmod_mat_clear(AT);
        nmod_mat_clear(AT2);
        nmod_mat_clear(B);
        nmod_mat_clear(BT);
    }

    /* Windows: row strides larger than the number of columns, and nothing
       written outside the destination window */
    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A, B, Aw, Bw, B2;
        slong r0, c0, r1, c1, i, j;

        m = rand_dim(state);
        n = rand_dim(state);
        r0 = n_randint(state, 3);
        c0 = n_randint(state, 9);
        r1 = n_randint(state, 3);
        c1 = n_randint(state, 9);

        mod = n_randtest_not_zero(state);

        nmod_mat_init(A, m + r0 + 2, n + c0 + 9, mod);
        nmod_mat_init(B, n + r1 + 2, m + c1 + 9, mod);
        nmod_mat_init(B2, n + r1 + 2, m + c1 + 9, mod);

        nmod_mat_randtest(A, state);
        nmod_mat_randtest(B, state);
        nmod_mat_set(B2, B);

        nmod_mat_window_init(Aw, A, r0, c0, r0 + m, c0 + n);
        nmod_mat_window_init(Bw, B, r1, c1, r1 + n, c1 + m);

        nmod_mat_transpose(Bw, Aw);

        if (!is_transpose(Bw, Aw))
            TEST_FUNCTION_FAIL("window: Bw != Aw^T\n");

        /* the same, entrywise, on a copy of B */
        for (i = 0; i < n; i++)
            for (j = 0; j < m; j++)
                nmod_mat_entry(B2, r1 + i, c1 + j) = nmod_mat_entry(Aw, j, i);

        if (!nmod_mat_equal(B, B2))
            TEST_FUNCTION_FAIL("window: write outside the destination\n");

        nmod_mat_window_clear(Aw);
        nmod_mat_window_clear(Bw);
        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(B2);
    }

    /* _nmod_mat_transpose on plain arrays, with padding after each row that
       must be left untouched */
    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        nn_ptr a, b, b2;
        slong Astride, Bstride, i, j;

        m = rand_dim(state);
        n = rand_dim(state);
        Astride = n + n_randint(state, 9);
        Bstride = m + n_randint(state, 9);

        a = flint_malloc((m * Astride + 1) * sizeof(ulong));
        b = flint_malloc((n * Bstride + 1) * sizeof(ulong));
        b2 = flint_malloc((n * Bstride + 1) * sizeof(ulong));

        for (i = 0; i < m * Astride; i++)
            a[i] = n_randtest(state);
        for (i = 0; i < n * Bstride; i++)
            b[i] = b2[i] = n_randtest(state);

        _nmod_mat_transpose(b, Bstride, a, Astride, m, n);

        for (i = 0; i < n; i++)
            for (j = 0; j < m; j++)
                b2[i * Bstride + j] = a[j * Astride + i];

        for (i = 0; i < n * Bstride; i++)
            if (b[i] != b2[i])
                TEST_FUNCTION_FAIL("_nmod_mat_transpose: m = %wd, n = %wd, "
                    "Astride = %wd, Bstride = %wd, index %wd\n",
                    m, n, Astride, Bstride, i);

        flint_free(a);
        flint_free(b);
        flint_free(b2);
    }

    /* Self-transpose, also of a window */
    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A, B, Bw, C;
        slong r0, c0;

        m = rand_dim(state);
        r0 = n_randint(state, 3);
        c0 = n_randint(state, 9);
        mod = n_randtest_not_zero(state);

        nmod_mat_init(A, m, m, mod);
        nmod_mat_init(B, m + r0 + 2, m + c0 + 9, mod);
        nmod_mat_init(C, m + r0 + 2, m + c0 + 9, mod);

        nmod_mat_randtest(B, state);
        nmod_mat_set(C, B);
        nmod_mat_window_init(Bw, B, r0, c0, r0 + m, c0 + m);
        nmod_mat_set(A, Bw);

        nmod_mat_transpose(Bw, Bw);

        if (!is_transpose(Bw, A))
            TEST_FUNCTION_FAIL("B != A^T\n");

        nmod_mat_transpose(Bw, Bw);

        if (!nmod_mat_equal(B, C))
            TEST_FUNCTION_FAIL("B != A\n");

        nmod_mat_window_clear(Bw);
        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
    }

    TEST_FUNCTION_END(state);
}
