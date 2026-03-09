/*
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_snf_transform, state)
{
    slong iter;

    /* Randomized tests: square, non-square, rank-deficient */
    for (iter = 0; iter < 3000 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A, S, S2, U, V, T1, T2;
        fmpz_t det;
        slong m, n, r, b, d, i;

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        r = n_randint(state, FLINT_MIN(m, n) + 1);
        d = FLINT_MIN(m, n);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(S, m, n);
        fmpz_mat_init(S2, m, n);
        fmpz_mat_init(U, m, m);
        fmpz_mat_init(V, n, n);
        fmpz_mat_init(T1, m, n);
        fmpz_mat_init(T2, m, n);
        fmpz_init(det);

        b = 1 + n_randint(state, 10) * n_randint(state, 10);
        fmpz_mat_randrank(A, state, r, b);

        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, n_randint(state, 2 * m * n + 1));

        fmpz_mat_snf_transform(S, U, V, A);

        /* Check S is in SNF */
        if (!fmpz_mat_is_in_snf(S))
        {
            flint_printf("FAIL:\n");
            flint_printf("S not in SNF, m=%wd n=%wd r=%wd\n", m, n, r);
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(S); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        /* Check U * A * V == S */
        fmpz_mat_mul(T1, U, A);
        fmpz_mat_mul(T2, T1, V);
        if (!fmpz_mat_equal(T2, S))
        {
            flint_printf("FAIL:\n");
            flint_printf("U*A*V != S, m=%wd n=%wd r=%wd\n", m, n, r);
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(S); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        /* Check |det(U)| == 1 (unimodular) */
        if (m > 0)
        {
            fmpz_mat_det(det, U);
            fmpz_abs(det, det);
            if (!fmpz_is_one(det))
            {
                flint_printf("FAIL:\n");
                flint_printf("|det(U)| != 1, m=%wd n=%wd\n", m, n);
                fflush(stdout);
                flint_abort();
            }
        }

        /* Check |det(V)| == 1 (unimodular) */
        if (n > 0)
        {
            fmpz_mat_det(det, V);
            fmpz_abs(det, det);
            if (!fmpz_is_one(det))
            {
                flint_printf("FAIL:\n");
                flint_printf("|det(V)| != 1, m=%wd n=%wd\n", m, n);
                fflush(stdout);
                flint_abort();
            }
        }

        /* Check rank: number of nonzero diagonal entries */
        {
            slong snf_rank = 0;
            for (i = 0; i < d; i++)
                if (!fmpz_is_zero(fmpz_mat_entry(S, i, i)))
                    snf_rank++;

            if (snf_rank != fmpz_mat_rank(A))
            {
                flint_printf("FAIL:\n");
                flint_printf("snf rank %wd != matrix rank %wd\n",
                    snf_rank, fmpz_mat_rank(A));
                fflush(stdout);
                flint_abort();
            }
        }

        /* Cross-check against fmpz_mat_snf */
        fmpz_mat_snf(S2, A);
        if (!fmpz_mat_equal(S, S2))
        {
            flint_printf("FAIL:\n");
            flint_printf("snf_transform disagrees with snf\n");
            flint_printf("m=%wd n=%wd r=%wd\n", m, n, r);
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            flint_printf("snf_transform: ");
            fmpz_mat_print_pretty(S); flint_printf("\n\n");
            flint_printf("snf:           ");
            fmpz_mat_print_pretty(S2); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(det);
        fmpz_mat_clear(T2);
        fmpz_mat_clear(T1);
        fmpz_mat_clear(V);
        fmpz_mat_clear(U);
        fmpz_mat_clear(S2);
        fmpz_mat_clear(S);
        fmpz_mat_clear(A);
    }

    /* Edge cases */
    {
        /* 0x0 */
        fmpz_mat_t S, U, V;
        fmpz_mat_init(S, 0, 0);
        fmpz_mat_init(U, 0, 0);
        fmpz_mat_init(V, 0, 0);
        fmpz_mat_snf_transform(S, U, V, S);
        fmpz_mat_clear(V);
        fmpz_mat_clear(U);
        fmpz_mat_clear(S);
    }
    {
        /* 0x5 */
        fmpz_mat_t A, S, U, V;
        fmpz_mat_init(A, 0, 5);
        fmpz_mat_init(S, 0, 5);
        fmpz_mat_init(U, 0, 0);
        fmpz_mat_init(V, 5, 5);
        fmpz_mat_snf_transform(S, U, V, A);
        fmpz_mat_clear(V);
        fmpz_mat_clear(U);
        fmpz_mat_clear(S);
        fmpz_mat_clear(A);
    }
    {
        /* 3x0 */
        fmpz_mat_t A, S, U, V;
        fmpz_mat_init(A, 3, 0);
        fmpz_mat_init(S, 3, 0);
        fmpz_mat_init(U, 3, 3);
        fmpz_mat_init(V, 0, 0);
        fmpz_mat_snf_transform(S, U, V, A);
        fmpz_mat_clear(V);
        fmpz_mat_clear(U);
        fmpz_mat_clear(S);
        fmpz_mat_clear(A);
    }
    {
        /* 1x1 negative */
        fmpz_mat_t A, S, U, V, T1, T2;
        fmpz_mat_init(A, 1, 1);
        fmpz_mat_init(S, 1, 1);
        fmpz_mat_init(U, 1, 1);
        fmpz_mat_init(V, 1, 1);
        fmpz_mat_init(T1, 1, 1);
        fmpz_mat_init(T2, 1, 1);

        fmpz_set_si(fmpz_mat_entry(A, 0, 0), -7);
        fmpz_mat_snf_transform(S, U, V, A);

        if (!fmpz_equal_si(fmpz_mat_entry(S, 0, 0), 7))
        {
            flint_printf("FAIL: 1x1 negative\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_mul(T1, U, A);
        fmpz_mat_mul(T2, T1, V);
        if (!fmpz_mat_equal(T2, S))
        {
            flint_printf("FAIL: U*A*V != S for 1x1\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(T2);
        fmpz_mat_clear(T1);
        fmpz_mat_clear(V);
        fmpz_mat_clear(U);
        fmpz_mat_clear(S);
        fmpz_mat_clear(A);
    }
    {
        /* Identity matrix */
        fmpz_mat_t A, S, U, V;
        fmpz_mat_init(A, 5, 5);
        fmpz_mat_init(S, 5, 5);
        fmpz_mat_init(U, 5, 5);
        fmpz_mat_init(V, 5, 5);

        fmpz_mat_one(A);
        fmpz_mat_snf_transform(S, U, V, A);

        if (!fmpz_mat_is_one(S))
        {
            flint_printf("FAIL: identity\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(V);
        fmpz_mat_clear(U);
        fmpz_mat_clear(S);
        fmpz_mat_clear(A);
    }
    {
        /* Zero matrix */
        fmpz_mat_t A, S, U, V;
        fmpz_mat_init(A, 4, 6);
        fmpz_mat_init(S, 4, 6);
        fmpz_mat_init(U, 4, 4);
        fmpz_mat_init(V, 6, 6);

        fmpz_mat_zero(A);
        fmpz_mat_snf_transform(S, U, V, A);

        if (!fmpz_mat_is_zero(S))
        {
            flint_printf("FAIL: zero matrix\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(V);
        fmpz_mat_clear(U);
        fmpz_mat_clear(S);
        fmpz_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}
