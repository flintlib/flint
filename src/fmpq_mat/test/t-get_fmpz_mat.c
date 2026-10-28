/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpq.h"
#include "fmpq_mat.h"

/* rowwise, colwise and matwise clearing of denominators, on windows,
   against direct computation */
TEST_FUNCTION_START(fmpq_mat_get_fmpz_mat, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        fmpq_mat_t A, W;
        fmpz_mat_t N, N2, M, M2;
        fmpz * den, * den2;
        fmpz_t D, d, t;
        slong r, c, i, j, r1, c1, bc;
        int shape, which;

        r = n_randint(state, 8);
        c = n_randint(state, 8);
        r1 = n_randint(state, 3);
        c1 = n_randint(state, 3);
        shape = n_randint(state, 3);
        which = n_randint(state, 4);
        bc = n_randint(state, 3);

        fmpq_mat_init(A, r + r1 + 1, c + c1 + 1);
        fmpq_mat_randtest(A, state, 1 + n_randint(state, 150));
        fmpz_init(D);
        fmpz_init(d);
        fmpz_init(t);
        fmpz_randtest_not_zero(D, state, 1 + n_randint(state, 150));
        fmpz_abs(D, D);

        /* structured denominators */
        if (shape != 0)
        {
            for (i = 0; i < A->r; i++)
            {
                for (j = 0; j < A->c; j++)
                {
                    if (shape == 1)
                        fmpz_set(fmpq_mat_entry_den(A, i, j), n_randint(state, 3) ? D : t);
                    else
                        fmpz_set_ui(fmpq_mat_entry_den(A, i, j), 1 + n_randint(state, 12));
                    if (fmpz_is_zero(fmpq_mat_entry_den(A, i, j)))
                        fmpz_one(fmpq_mat_entry_den(A, i, j));
                    fmpq_canonicalise(fmpq_mat_entry(A, i, j));
                }
            }
        }

        fmpq_mat_window_init(W, A, r1, c1, r1 + r, c1 + c);
        fmpz_mat_init(N, r, c);
        fmpz_mat_init(N2, r, c);
        fmpz_mat_init(M, r, bc);
        fmpz_mat_init(M2, r, bc);
        den = _fmpz_vec_init(FLINT_MAX(r, c) + 1);
        den2 = _fmpz_vec_init(FLINT_MAX(r, c) + 1);

        if (which == 0)
        {
            fmpq_mat_get_fmpz_mat_rowwise(N, den, W);
            for (i = 0; i < r && c > 0; i++)
            {
                fmpz_one(den2 + i);
                for (j = 0; j < c; j++)
                    fmpz_lcm(den2 + i, den2 + i, fmpq_mat_entry_den(W, i, j));
                for (j = 0; j < c; j++)
                {
                    fmpz_divexact(t, den2 + i, fmpq_mat_entry_den(W, i, j));
                    fmpz_mul(fmpz_mat_entry(N2, i, j), fmpq_mat_entry_num(W, i, j), t);
                }
            }
            if (c > 0 && (!_fmpz_vec_equal(den, den2, r) || !fmpz_mat_equal(N, N2)))
                TEST_FUNCTION_FAIL("rowwise: r = %wd, c = %wd, shape = %d\n", r, c, shape);
        }
        else if (which == 1)
        {
            fmpq_mat_get_fmpz_mat_colwise(N, den, W);
            for (j = 0; j < c && r > 0; j++)
            {
                fmpz_one(den2 + j);
                for (i = 0; i < r; i++)
                    fmpz_lcm(den2 + j, den2 + j, fmpq_mat_entry_den(W, i, j));
                for (i = 0; i < r; i++)
                {
                    fmpz_divexact(t, den2 + j, fmpq_mat_entry_den(W, i, j));
                    fmpz_mul(fmpz_mat_entry(N2, i, j), fmpq_mat_entry_num(W, i, j), t);
                }
            }
            if (r > 0 && (!_fmpz_vec_equal(den, den2, c) || !fmpz_mat_equal(N, N2)))
                TEST_FUNCTION_FAIL("colwise: r = %wd, c = %wd, shape = %d\n", r, c, shape);
        }
        else if (which == 2)
        {
            fmpq_mat_get_fmpz_mat_matwise(N, d, W);
            fmpz_one(den2);
            for (i = 0; i < r; i++)
                for (j = 0; j < c; j++)
                    fmpz_lcm(den2, den2, fmpq_mat_entry_den(W, i, j));
            for (i = 0; i < r; i++)
                for (j = 0; j < c; j++)
                {
                    fmpz_divexact(t, den2, fmpq_mat_entry_den(W, i, j));
                    fmpz_mul(fmpz_mat_entry(N2, i, j), fmpq_mat_entry_num(W, i, j), t);
                }
            if (!fmpz_equal(d, den2) || !fmpz_mat_equal(N, N2))
                TEST_FUNCTION_FAIL("matwise: r = %wd, c = %wd, shape = %d\n", r, c, shape);
        }
        else
        {
            /* two matrices, the second of 0 to 2 columns */
            fmpq_mat_t B;
            fmpq_mat_init(B, r, bc);
            fmpq_mat_randtest(B, state, 50);
            fmpq_mat_get_fmpz_mat_rowwise_2(N, M, den, W, B);
            for (i = 0; i < r && c > 0; i++)
            {
                fmpz_one(den2 + i);
                for (j = 0; j < c; j++)
                    fmpz_lcm(den2 + i, den2 + i, fmpq_mat_entry_den(W, i, j));
                for (j = 0; j < bc; j++)
                    fmpz_lcm(den2 + i, den2 + i, fmpq_mat_entry_den(B, i, j));
                for (j = 0; j < c; j++)
                {
                    fmpz_divexact(t, den2 + i, fmpq_mat_entry_den(W, i, j));
                    fmpz_mul(fmpz_mat_entry(N2, i, j), fmpq_mat_entry_num(W, i, j), t);
                }
                for (j = 0; j < bc; j++)
                {
                    fmpz_divexact(t, den2 + i, fmpq_mat_entry_den(B, i, j));
                    fmpz_mul(fmpz_mat_entry(M2, i, j), fmpq_mat_entry_num(B, i, j), t);
                }
            }
            if (c > 0 && (!_fmpz_vec_equal(den, den2, r) || !fmpz_mat_equal(N, N2) || !fmpz_mat_equal(M, M2)))
                TEST_FUNCTION_FAIL("rowwise_2: r = %wd, c = %wd, shape = %d\n", r, c, shape);
            fmpq_mat_clear(B);
        }

        fmpq_mat_window_clear(W);
        fmpq_mat_clear(A);
        fmpz_mat_clear(N);
        fmpz_mat_clear(N2);
        fmpz_mat_clear(M);
        fmpz_mat_clear(M2);
        _fmpz_vec_clear(den, FLINT_MAX(r, c) + 1);
        _fmpz_vec_clear(den2, FLINT_MAX(r, c) + 1);
        fmpz_clear(D);
        fmpz_clear(d);
        fmpz_clear(t);
    }

    TEST_FUNCTION_END(state);
}
