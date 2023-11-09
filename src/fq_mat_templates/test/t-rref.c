/*
    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"
#include "perm.h"

int
check_rref_form(slong * perm, TEMPLATE(T, mat_t) A, slong rank,
                const TEMPLATE(T, ctx_t) ctx)
{
    slong i, j, k, prev_pivot;

    /* bottom should be zero */
    for (i = rank; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            if (!TEMPLATE(T, is_zero) (TEMPLATE(T, mat_entry) (A, i, j), ctx))
                return 0;

    prev_pivot = -1;

    for (i = 0; i < rank; i++)
    {
        for (j = 0; j < A->c; j++)
        {
            if (!TEMPLATE(T, is_zero) (TEMPLATE(T, mat_entry) (A, i, j), ctx))
            {
                /* pivot should have a higher column index than previous */
                if (j <= prev_pivot)
                    return 0;

                /* column should be 0 ... 0 1 0 ... 0 */
                for (k = 0; k < rank; k++)
                {
                    if (i == k)
                    {
                        if (!TEMPLATE(T, is_one)
                            (TEMPLATE(T, mat_entry) (A, k, j), ctx))
                            return 0;
                    }
                    else
                    {
                        if (!TEMPLATE(T, is_zero)
                            (TEMPLATE(T, mat_entry) (A, k, j), ctx))
                            return 0;
                    }
                }

                prev_pivot = j;
                break;
            }
        }
    }

    return 1;
}

TEST_TEMPLATE_FUNCTION_START(T, mat_rref, state)
{
    slong i;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, mat_t) A, B, C, D;
        TEMPLATE(T, t) c;
        slong j, k, m, n, rank1, rank2;
        slong *perm;
        int equal;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, init) (c, ctx);

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        perm = _perm_init(2 * m);

        TEMPLATE(T, mat_init) (A, m, n, ctx);
        TEMPLATE(T, mat_init) (D, 2 * m, n, ctx);

        TEMPLATE(T, mat_randtest) (A, state, ctx);
        TEMPLATE(T, mat_init_set) (B, A, ctx);
        TEMPLATE(T, mat_init_set) (C, A, ctx);

        rank1 = TEMPLATE(T, mat_rref) (B, ctx);

        if (!check_rref_form(perm, B, rank1, ctx))
        {
            printf("FAIL (malformed rref)\n");
            TEMPLATE(T, mat_print_pretty) (A, ctx);
            printf("\n\n");
            TEMPLATE(T, mat_print_pretty) (B, ctx);
            printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        /* Concatenate the original matrix with the rref, scramble the rows,
           and check that the rref is the same */
        _perm_randtest(perm, 2 * m, state);

        for (j = 0; j < m; j++)
        {
            TEMPLATE(T, randtest_not_zero) (c, state, ctx);
            for (k = 0; k < n; k++)
            {
                TEMPLATE(T, mul) (TEMPLATE(T, mat_entry) (D, perm[j], k),
                                  TEMPLATE(T, mat_entry) (A, j, k), c, ctx);
            }
        }

        for (j = 0; j < m; j++)
        {
            TEMPLATE(T, randtest_not_zero) (c, state, ctx);
            for (k = 0; k < n; k++)
            {
                TEMPLATE(T, mul) (TEMPLATE(T, mat_entry) (D, perm[m + j], k),
                                  TEMPLATE(T, mat_entry) (B, j, k), c, ctx);

            }
        }

        rank2 = TEMPLATE(T, mat_rref) (D, ctx);
        equal = (rank1 == rank2);

        if (equal)
        {
            for (j = 0; j < rank2; j++)
                for (k = 0; k < n; k++)
                {
                    equal = equal
                        && TEMPLATE(T,
                                    equal) (TEMPLATE(T, mat_entry) (B, j, k),
                                            TEMPLATE(T, mat_entry) (D, j, k),
                                            ctx);
                }
            for (j = rank2; j < 2 * rank2; j++)
                for (k = 0; k < n; k++)
                {
                    equal = equal
                        && TEMPLATE(T,
                                    is_zero) (TEMPLATE(T, mat_entry) (D, j, k),
                                              ctx);
                }
        }

        if (!equal)
        {
            flint_printf("FAIL (rank1 = %wd, rank2 = %wd)!\n", rank1, rank2);
            TEMPLATE(T, mat_print_pretty) (A, ctx);
            printf("\n\n");
            TEMPLATE(T, mat_print_pretty) (B, ctx);
            printf("\n\n");
            TEMPLATE(T, mat_print_pretty) (D, ctx);
            printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        _perm_clear(perm);
        TEMPLATE(T, mat_clear) (A, ctx);
        TEMPLATE(T, mat_clear) (B, ctx);
        TEMPLATE(T, mat_clear) (C, ctx);
        TEMPLATE(T, mat_clear) (D, ctx);

        TEMPLATE(T, clear) (c, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
