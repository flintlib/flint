/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod.h"
#include "nmod_mat.h"

/* Defined in t-mul.c and t-mul_classical_threaded.c */
#ifndef nmod_mat_mul_check
#define nmod_mat_mul_check nmod_mat_mul_check
void
nmod_mat_mul_check(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    slong i, j, k;

    ulong s0, s1, s2;
    ulong t0, t1;

    for (i = 0; i < A->r; i++)
    {
        for (j = 0; j < B->c; j++)
        {
            s0 = s1 = s2 = UWORD(0);

            for (k = 0; k < A->c; k++)
            {
                umul_ppmm(t1, t0, nmod_mat_entry(A, i, k), nmod_mat_entry(B, k, j));
                add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, t1, t0);
            }

            NMOD_RED(s2, s2, C->mod);
            NMOD_RED3(s0, s2, s1, s0, C->mod);
            nmod_mat_entry(C, i, j) = s0;
        }
    }
}
#endif

TEST_FUNCTION_START(nmod_mat_mul, state)
{
    slong i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, B, C, D;
        ulong mod;

        slong m, k, n;

        m = n_randint(state, 75);
        k = n_randint(state, 75);
        n = n_randint(state, 75);

        /* We want to generate matrices with many entries close to half
           or full limbs with high probability, to stress overflow handling */
        switch (n_randint(state, 3))
        {
            case 0:
                mod = n_randtest_not_zero(state);
                break;
            case 1:
                mod = UWORD_MAX/2 + 1 - n_randbits(state, 4);
                break;
            case 2:
            default:
                mod = UWORD_MAX - n_randbits(state, 4);
                break;
        }

        nmod_mat_init(A, m, n, mod);
        nmod_mat_init(B, n, k, mod);
        nmod_mat_init(C, m, k, mod);
        nmod_mat_init(D, m, k, mod);

        if (n_randint(state, 2))
            nmod_mat_randtest(A, state);
        else
            nmod_mat_randfull(A, state);

        if (n_randint(state, 2))
            nmod_mat_randtest(B, state);
        else
            nmod_mat_randfull(B, state);

        nmod_mat_randtest(C, state);  /* make sure noise in the output is ok */

        nmod_mat_mul(C, A, B);
        nmod_mat_mul_check(D, A, B);

        if (!nmod_mat_equal(C, D))
            TEST_FUNCTION_FAIL(
                    "Results not equal\n"
                    "A = %{nmod_mat}\n"
                    "B = %{nmod_mat}\n"
                    "C = %{nmod_mat}\n"
                    "D = %{nmod_mat}\n",
                    A, B, C, D);

        if (n == k)
        {
            nmod_mat_mul(A, A, B);

            if (!nmod_mat_equal(A, C))
                TEST_FUNCTION_FAIL("Aliasing failed\n");
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        nmod_mat_clear(D);
    }

    /* thin shapes (one dimension up to 8, the others larger), which
       nmod_mat_mul sends to specific code, on windows of larger matrices
       half of the time */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, B, C, D, PA, PB, PC, PC0;
        slong m, k, n, ra, ca, rb, cb, rc, cc, r, c;
        ulong mod;
        int win = n_randint(state, 2);

        m = 1 + n_randint(state, 150);
        k = 1 + n_randint(state, 150);
        n = 1 + n_randint(state, 150);
        switch (n_randint(state, 4))
        {
            case 0: m = 1 + n_randint(state, 8); break;
            case 1: k = 1 + n_randint(state, 8); break;
            case 2: n = 1 + n_randint(state, 8); break;
            default: n = 1; break;
        }

        if (n_randint(state, 4) == 0)
            mod = n_randtest_not_zero(state);
        else
            mod = FLINT_MAX(UWORD(1), n_randbits(state, 1 + n_randint(state, FLINT_MIN(52, FLINT_BITS))));

        ra = win ? n_randint(state, 4) : 0;
        ca = win ? n_randint(state, 4) : 0;
        rb = win ? n_randint(state, 4) : 0;
        cb = win ? n_randint(state, 4) : 0;
        rc = win ? n_randint(state, 4) : 0;
        cc = win ? n_randint(state, 4) : 0;

        nmod_mat_init(PA, ra + m + (win ? n_randint(state, 4) : 0),
                          ca + k + (win ? n_randint(state, 4) : 0), mod);
        nmod_mat_init(PB, rb + k + (win ? n_randint(state, 4) : 0),
                          cb + n + (win ? n_randint(state, 4) : 0), mod);
        nmod_mat_init(PC, rc + m + (win ? n_randint(state, 4) : 0),
                          cc + n + (win ? n_randint(state, 4) : 0), mod);
        nmod_mat_randfull(PA, state);
        nmod_mat_randfull(PB, state);
        nmod_mat_randtest(PC, state);
        nmod_mat_init_set(PC0, PC);

        nmod_mat_window_init(A, PA, ra, ca, ra + m, ca + k);
        nmod_mat_window_init(B, PB, rb, cb, rb + k, cb + n);
        nmod_mat_window_init(C, PC, rc, cc, rc + m, cc + n);
        nmod_mat_init(D, m, n, mod);

        flint_set_num_threads(n_randint(state, 2) ? 1 : 1 + n_randint(state, 4));

        nmod_mat_mul(C, A, B);
        nmod_mat_mul_check(D, A, B);

        if (!nmod_mat_equal(C, D))
            TEST_FUNCTION_FAIL("thin shapes: m: %wd, k: %wd, n: %wd, "
                               "mod: %wu, windows: %d, threads: %d\n",
                               m, k, n, mod, win, flint_get_num_threads());

        for (r = 0; r < PC->r; r++)
            for (c = 0; c < PC->c; c++)
                if ((r < rc || r >= rc + m || c < cc || c >= cc + n)
                        && nmod_mat_entry(PC, r, c) != nmod_mat_entry(PC0, r, c))
                    TEST_FUNCTION_FAIL("thin shapes: write outside the "
                                       "window: m: %wd, k: %wd, n: %wd, "
                                       "mod: %wu, threads: %d\n", m, k, n,
                                       mod, flint_get_num_threads());

        /* aliasing, when the shapes allow it */
        if (k == n)
        {
            nmod_mat_mul(A, A, B);
            if (!nmod_mat_equal(A, D))
                TEST_FUNCTION_FAIL("thin shapes: aliasing C = A: m: %wd, "
                                   "k: %wd, mod: %wu\n", m, k, mod);
        }
        else if (m == k)
        {
            nmod_mat_mul(B, A, B);
            if (!nmod_mat_equal(B, D))
                TEST_FUNCTION_FAIL("thin shapes: aliasing C = B: m: %wd, "
                                   "n: %wd, mod: %wu\n", m, n, mod);
        }

        nmod_mat_window_clear(A);
        nmod_mat_window_clear(B);
        nmod_mat_window_clear(C);
        nmod_mat_clear(D);
        nmod_mat_clear(PA);
        nmod_mat_clear(PB);
        nmod_mat_clear(PC);
        nmod_mat_clear(PC0);
    }

    /* Test aliasing with windows */
    {
        nmod_mat_t A, B, A_window;

        nmod_mat_init(A, 2, 2, 3);
        nmod_mat_init(B, 2, 2, 3);

        nmod_mat_window_init(A_window, A, 0, 0, 2, 2);

        nmod_mat_one(A);
        nmod_mat_one(B);
        nmod_mat_entry(B, 0, 1) = 1;
        nmod_mat_entry(B, 1, 0) = 1;

        nmod_mat_mul(A_window, B, A_window);

        if (!nmod_mat_equal(A, B))
            TEST_FUNCTION_FAIL(
                    "Window aliasing failed\n"
                    "A = %{nmod_mat}\n"
                    "B = %{nmod_mat}\n",
                    A, B);

        nmod_mat_window_clear(A_window);
        nmod_mat_clear(A);
        nmod_mat_clear(B);
    }

    TEST_FUNCTION_END(state);
}
