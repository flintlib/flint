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
#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_mat.h"

/* checks D = A0 + x B0 entry by entry (A0 and B0: copies of the inputs) */
static int
_scalar_addmul_check(const nmod_mat_t D, const nmod_mat_t A0,
                     const nmod_mat_t B0, ulong x)
{
    slong i, j;
    ulong c = nmod_set_ui(x, D->mod);

    for (i = 0; i < D->r; i++)
        for (j = 0; j < D->c; j++)
            if (nmod_mat_entry(D, i, j) != nmod_add(nmod_mat_entry(A0, i, j),
                    nmod_mul(c, nmod_mat_entry(B0, i, j), D->mod), D->mod))
                return 0;
    return 1;
}

TEST_FUNCTION_START(nmod_mat_scalar_addmul_ui, state)
{
    slong m, n, rep;
    ulong mod, x;

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A, B, D, A0, B0, big, W;
        slong r0, c0;
        int alias;

        m = n_randint(state, 40);
        n = n_randint(state, 40);
        mod = n_randtest_not_zero(state);
        /* any word, and the special values 0, 1, n - 1 */
        switch (n_randint(state, 5))
        {
            case 0: x = 0; break;
            case 1: x = 1; break;
            case 2: x = mod - 1; break;
            case 3: x = n_randlimb(state); break;
            default: x = n_randint(state, mod);
        }

        nmod_mat_init(A, m, n, mod);
        nmod_mat_init(B, m, n, mod);
        nmod_mat_init(D, m, n, mod);
        nmod_mat_randtest(A, state);
        nmod_mat_randtest(B, state);
        nmod_mat_init_set(A0, A);
        nmod_mat_init_set(B0, B);

        /* D = A + x B, D distinct from A and B, then aliased with A, B or
           both */
        alias = n_randint(state, 4);
        if (alias == 0)
        {
            nmod_mat_scalar_addmul_ui(D, A, B, x);
            if (!_scalar_addmul_check(D, A0, B0, x))
                TEST_FUNCTION_FAIL("A = %{nmod_mat}\nB = %{nmod_mat}\nx = %wu\nD = %{nmod_mat}\n",
                                   A0, B0, x, D);
        }
        else if (alias == 1)
        {
            nmod_mat_scalar_addmul_ui(A, A, B, x);
            if (!_scalar_addmul_check(A, A0, B0, x))
                TEST_FUNCTION_FAIL("aliasing C = A: x = %wu\n", x);
        }
        else if (alias == 2)
        {
            nmod_mat_scalar_addmul_ui(B, A, B, x);
            if (!_scalar_addmul_check(B, A0, B0, x))
                TEST_FUNCTION_FAIL("aliasing C = B: x = %wu\n", x);
        }
        else
        {
            nmod_mat_scalar_addmul_ui(A, A, A, x);
            if (!_scalar_addmul_check(A, A0, A0, x))
                TEST_FUNCTION_FAIL("aliasing C = A = B: x = %wu\n", x);
        }

        /* windows: rows that do not follow each other */
        nmod_mat_init(big, m + 2, n + 3, mod);
        nmod_mat_randtest(big, state);
        r0 = n_randint(state, 3);
        c0 = n_randint(state, 4);
        nmod_mat_window_init(W, big, r0, c0, r0 + m, c0 + n);
        nmod_mat_scalar_addmul_ui(W, A0, B0, x);
        if (!_scalar_addmul_check(W, A0, B0, x))
            TEST_FUNCTION_FAIL("window: x = %wu\n", x);
        nmod_mat_set(D, W);
        nmod_mat_scalar_addmul_ui(W, W, B0, x);
        if (!_scalar_addmul_check(W, D, B0, x))
            TEST_FUNCTION_FAIL("window, aliasing C = A: x = %wu\n", x);
        nmod_mat_window_clear(W);

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(D);
        nmod_mat_clear(A0);
        nmod_mat_clear(B0);
        nmod_mat_clear(big);
    }

    TEST_FUNCTION_END(state);
}
