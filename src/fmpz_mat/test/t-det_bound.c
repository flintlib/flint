/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_det_bound, state)
{
    fmpz_mat_t A;
    slong i, m;

    fmpz_t det, bound;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 10);

        fmpz_mat_init(A, m, m);

        fmpz_init(det);
        fmpz_init(bound);

        fmpz_mat_randtest(A, state, 1+n_randint(state,200));

        fmpz_mat_det(det, A);
        fmpz_mat_det_bound(bound, A);

        if (fmpz_cmp(det, bound) > 0)
        {
            flint_printf("FAIL:\n");
            flint_printf("bound too small!\n");
            fmpz_mat_print_pretty(A), flint_printf("\n");
            flint_printf("det: "), fmpz_print(det), flint_printf("\n");
            flint_printf("bound: "), fmpz_print(bound), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(det);
        fmpz_clear(bound);
        fmpz_mat_clear(A);
    }

    /* Test fmpz_mat_det_bound_nonzero on random rectangular matrices.
       For each test matrix, extract a few random square submatrices,
       compute their determinants, and verify all are bounded by the
       nonzero bound. Also verify the bound is always positive.       */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t B, sub;
        fmpz_t det, bound;
        slong r, c, k, s, t, trial;
        slong *row_idx, *col_idx;

        r = n_randint(state, 8) + 1;
        c = n_randint(state, 8) + 1;

        fmpz_mat_init(B, r, c);
        fmpz_init(det);
        fmpz_init(bound);

        /* Occasionally insert zero rows/columns to exercise the
           nonzero-exclusion logic specifically                    */
        if (n_randint(state, 3) == 0)
            fmpz_mat_randtest(B, state, 1 + n_randint(state, 100));
        else
            fmpz_mat_randtest(B, state, 1 + n_randint(state, 100));

        if (n_randint(state, 4) == 0)
        {
            /* Zero out a random row */
            slong zr = n_randint(state, r);
            for (t = 0; t < c; t++)
                fmpz_zero(fmpz_mat_entry(B, zr, t));
        }

        if (n_randint(state, 4) == 0)
        {
            /* Zero out a random column */
            slong zc = n_randint(state, c);
            for (s = 0; s < r; s++)
                fmpz_zero(fmpz_mat_entry(B, s, zc));
        }

        fmpz_mat_det_bound_submatrix(bound, B);

        /* Bound must always be positive */
        if (fmpz_sgn(bound) <= 0)
        {
            flint_printf("FAIL (det_bound_submatrix):\n");
            flint_printf("bound is not positive!\n");
            fmpz_mat_print_pretty(B), flint_printf("\n");
            flint_printf("bound: "), fmpz_print(bound), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        /* Extract and test a few random square submatrices of size
           k x k where k <= min(r, c)                               */
        k = FLINT_MIN(r, c);
        k = n_randint(state, k) + 1;

        row_idx = (slong *) flint_malloc(r * sizeof(slong));
        col_idx = (slong *) flint_malloc(c * sizeof(slong));

        for (trial = 0; trial < 3; trial++)
        {
            slong ri, ci;

            /* Choose k distinct rows at random via partial Fisher-Yates */
            for (s = 0; s < r; s++) row_idx[s] = s;
            for (s = 0; s < k; s++)
            {
                slong swap = s + n_randint(state, r - s);
                slong tmp = row_idx[s]; row_idx[s] = row_idx[swap]; row_idx[swap] = tmp;
            }

            /* Choose k distinct columns at random */
            for (s = 0; s < c; s++) col_idx[s] = s;
            for (s = 0; s < k; s++)
            {
                slong swap = s + n_randint(state, c - s);
                slong tmp = col_idx[s]; col_idx[s] = col_idx[swap]; col_idx[swap] = tmp;
            }

            fmpz_mat_init(sub, k, k);

            for (ri = 0; ri < k; ri++)
                for (ci = 0; ci < k; ci++)
                    fmpz_set(fmpz_mat_entry(sub, ri, ci),
                             fmpz_mat_entry(B, row_idx[ri], col_idx[ci]));

            fmpz_mat_det(det, sub);
            fmpz_abs(det, det);

            if (fmpz_cmp(det, bound) > 0)
            {
                flint_printf("FAIL (det_bound_submatrix):\n");
                flint_printf("bound too small for submatrix determinant!\n");
                flint_printf("A (%wd x %wd):\n", r, c);
                fmpz_mat_print_pretty(B), flint_printf("\n");
                flint_printf("submatrix size: %wd\n", k);
                flint_printf("|det(sub)|: "), fmpz_print(det), flint_printf("\n");
                flint_printf("bound: "), fmpz_print(bound), flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_mat_clear(sub);
        }

        flint_free(row_idx);
        flint_free(col_idx);

        fmpz_clear(det);
        fmpz_clear(bound);
        fmpz_mat_clear(B);
    }

    TEST_FUNCTION_END(state);
}
