/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod_mat.h"

static void
check_rref(fmpz_mod_mat_t A)
{
    slong i, j, prev_pivot, prev_row_zero;

    prev_row_zero = 0;
    prev_pivot = -1;

    for (i = 0; i < A->r; i++)
    {
        for (j = 0; j < A->c; j++)
        {
            /* Found nonzero entry */
            if (!fmpz_is_zero(fmpz_mod_mat_entry(A, i, j)))
            {
                if (prev_row_zero)
                {
                    flint_printf("nonzero row after zero row\n");
                    fflush(stdout);
                    flint_abort();
                }

                if (j <= prev_pivot)
                {
                    flint_printf("pivot not strictly to the right of previous\n");
                    fflush(stdout);
                    flint_abort();
                }

                prev_pivot = j;
                break;
            }

            prev_row_zero = (j + 1 == A->c);
        }
    }
}

TEST_FUNCTION_START(fmpz_mod_mat_rref, state)
{
    fmpz_mod_mat_t A;
    fmpz_mod_ctx_t ctx;
    slong i, m, n, d, r, rank;

    /* Maximally sparse matrices of given rank */
    for (i = 0; i < 10000; i++)
    {
        m = n_randint(state, 10);
        n = n_randint(state, 10);

        fmpz_mod_ctx_init_rand_bits_prime(ctx, state, 100);

        for (r = 0; r <= FLINT_MIN(m, n); r++)
        {
            fmpz_mod_mat_init(A, m, n, ctx);
            fmpz_mod_mat_randrank(A, state, r, ctx);

            rank = fmpz_mod_mat_rref(A, A, ctx);

            if (r < rank)
            {
                fmpz_mod_mat_print_pretty(A, ctx);
                flint_printf("FAIL:\n");
                flint_printf("wrong rank!\n");
                fflush(stdout);
                flint_abort();
            }

            check_rref(A);

            fmpz_mod_mat_clear(A, ctx);
        }

        fmpz_mod_ctx_clear(ctx);
    }

    /* Dense */
    for (i = 0; i < 10000; i++)
    {
        m = n_randint(state, 5);
        n = n_randint(state, 4);

        fmpz_mod_ctx_init_rand_bits_prime(ctx, state, 100);

        for (r = 0; r <= FLINT_MIN(m, n); r++)
        {
            d = n_randint(state, 2 * m * n + 1);

            fmpz_mod_mat_init(A, m, n, ctx);
            fmpz_mod_mat_randrank(A, state, r, ctx);
            fmpz_mod_mat_randops(A, state, d, ctx);

            rank = fmpz_mod_mat_rref(A, A, ctx);

            if (r < rank)
            {
                flint_printf("FAIL:\n");
                flint_printf("wrong rank!\n");
                fflush(stdout);
                flint_abort();
            }

            check_rref(A);

            fmpz_mod_mat_clear(A, ctx);
        }

        fmpz_mod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
