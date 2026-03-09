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
#include "fmpz_vec.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_elementary_divisors, state)
{
    slong iter;

    /* Randomized tests: cross-check against SNF diagonal */
    for (iter = 0; iter < 3000 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A, S;
        fmpz * ed;
        slong m, n, r, b, d, i, rank;

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        r = n_randint(state, FLINT_MIN(m, n) + 1);
        d = FLINT_MIN(m, n);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(S, m, n);

        b = 1 + n_randint(state, 10) * n_randint(state, 10);
        fmpz_mat_randrank(A, state, r, b);

        if (n_randint(state, 2))
            fmpz_mat_randops(A, state,
                n_randint(state, 2 * m * n + 1));

        ed = _fmpz_vec_init(d);
        fmpz_mat_elementary_divisors(ed, &rank, A);

        /* Check rank */
        if (rank != r)
        {
            flint_printf("FAIL:\n");
            flint_printf("rank %wd != expected %wd\n", rank, r);
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        /* Cross-check against SNF diagonal */
        fmpz_mat_snf(S, A);

        for (i = 0; i < r; i++)
        {
            if (!fmpz_equal(&ed[i], fmpz_mat_entry(S, i, i)))
            {
                flint_printf("FAIL:\n");
                flint_printf("ed[%wd] = %{fmpz} != snf[%wd] = %{fmpz}\n",
                    i, &ed[i], i, fmpz_mat_entry(S, i, i));
                flint_printf("m=%wd n=%wd r=%wd\n", m, n, r);
                fmpz_mat_print_pretty(A); flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }
        }

        _fmpz_vec_clear(ed, d);
        fmpz_mat_clear(S);
        fmpz_mat_clear(A);
    }

    /* Edge cases */
    {
        /* 0x0 matrix */
        fmpz_mat_t A;
        slong rank;
        fmpz_mat_init(A, 0, 0);
        fmpz_mat_elementary_divisors(NULL, &rank, A);
        if (rank != 0)
        {
            flint_printf("FAIL: 0x0 rank\n");
            fflush(stdout);
            flint_abort();
        }
        fmpz_mat_clear(A);
    }
    {
        /* Zero matrix */
        fmpz_mat_t A;
        slong rank;
        fmpz_mat_init(A, 5, 3);
        fmpz_mat_zero(A);
        fmpz_mat_elementary_divisors(NULL, &rank, A);
        if (rank != 0)
        {
            flint_printf("FAIL: zero matrix rank\n");
            fflush(stdout);
            flint_abort();
        }
        fmpz_mat_clear(A);
    }
    {
        /* Identity matrix */
        fmpz_mat_t A;
        fmpz * ed;
        slong rank, i;
        fmpz_mat_init(A, 4, 4);
        fmpz_mat_one(A);
        ed = _fmpz_vec_init(4);
        fmpz_mat_elementary_divisors(ed, &rank, A);
        if (rank != 4)
        {
            flint_printf("FAIL: identity rank\n");
            fflush(stdout);
            flint_abort();
        }
        for (i = 0; i < 4; i++)
        {
            if (!fmpz_is_one(&ed[i]))
            {
                flint_printf("FAIL: identity ed[%wd] = %{fmpz}\n",
                    i, &ed[i]);
                fflush(stdout);
                flint_abort();
            }
        }
        _fmpz_vec_clear(ed, 4);
        fmpz_mat_clear(A);
    }
    {
        /* Known diagonal matrix */
        fmpz_mat_t A;
        fmpz * ed;
        slong rank;
        fmpz_mat_init(A, 3, 3);
        fmpz_mat_zero(A);
        fmpz_set_si(fmpz_mat_entry(A, 0, 0), 12);
        fmpz_set_si(fmpz_mat_entry(A, 1, 1), 6);
        fmpz_set_si(fmpz_mat_entry(A, 2, 2), 0);
        ed = _fmpz_vec_init(3);
        fmpz_mat_elementary_divisors(ed, &rank, A);
        if (rank != 2)
        {
            flint_printf("FAIL: diag rank\n");
            fflush(stdout);
            flint_abort();
        }
        if (!fmpz_equal_si(&ed[0], 6)
            || !fmpz_equal_si(&ed[1], 12))
        {
            flint_printf("FAIL: diag ed = [%{fmpz}, %{fmpz}], "
                "expected [6, 12]\n", &ed[0], &ed[1]);
            fflush(stdout);
            flint_abort();
        }
        _fmpz_vec_clear(ed, 3);
        fmpz_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}
