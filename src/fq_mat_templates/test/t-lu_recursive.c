/*
    Copyright (C) 2011 Fredrik Johansson
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

/* Defined in t-lu_classical.c and t-lu_recursive.c */
#ifndef perm
#define perm perm
void
perm(TEMPLATE(T, mat_t) A, slong * P)
{
    slong i;
    TEMPLATE(T, struct) ** tmp;

    if (A->c == 0 || A->r == 0)
        return;

    tmp = flint_malloc(sizeof(TEMPLATE(T, struct) *) * A->r);

    for (i = 0; i < A->r; i++)
        tmp[P[i]] = A->rows[i];
    for (i = 0; i < A->r; i++)
        A->rows[i] = tmp[i];

    flint_free(tmp);
}
#endif

/* Defined in t-lu_classical.c and t-lu_recursive.c */
#ifndef check
#define check check
void
check(slong * P, TEMPLATE(T, mat_t) LU, const TEMPLATE(T, mat_t) A, slong rank,
      const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, mat_t) B, L, U;
    slong m, n, i, j;

    m = A->r;
    n = A->c;

    TEMPLATE(T, mat_init) (B, m, n, ctx);
    TEMPLATE(T, mat_init) (L, m, m, ctx);
    TEMPLATE(T, mat_init) (U, m, n, ctx);

    rank = FLINT_ABS(rank);

    for (i = rank; i < FLINT_MIN(m, n); i++)
    {
        for (j = i; j < n; j++)
        {
            if (!TEMPLATE(T, is_zero) (TEMPLATE(T, mat_entry) (LU, i, j), ctx))
            {
                printf("FAIL: wrong shape!\n");
                fflush(stdout);
                flint_abort();
            }
        }
    }

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < FLINT_MIN(i, n); j++)
            TEMPLATE(T, mat_entry_set) (L, i, j,
                                        TEMPLATE(T, mat_entry) (LU, i, j),
                                        ctx);
        if (i < rank)
            TEMPLATE(T, one) (TEMPLATE(T, mat_entry) (L, i, i), ctx);
        for (j = i; j < n; j++)
            TEMPLATE(T, mat_entry_set) (U, i, j,
                                        TEMPLATE(T, mat_entry) (LU, i, j),
                                        ctx);
    }

    TEMPLATE(T, mat_mul) (B, L, U, ctx);
    perm(B, P);

    if (!TEMPLATE(T, mat_equal) (A, B, ctx))
    {
        printf("FAIL\n");
        printf("A:\n");
        TEMPLATE(T, mat_print_pretty) (A, ctx);
        printf("LU:\n");
        TEMPLATE(T, mat_print_pretty) (LU, ctx);
        printf("B:\n");
        TEMPLATE(T, mat_print_pretty) (B, ctx);
        fflush(stdout);
        flint_abort();
    }

    TEMPLATE(T, mat_clear) (B, ctx);
    TEMPLATE(T, mat_clear) (L, ctx);
    TEMPLATE(T, mat_clear) (U, ctx);
}
#endif

TEST_TEMPLATE_FUNCTION_START(T, mat_lu_recursive, state)
{
    slong i;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, mat_t) A, LU;

        slong m, n, r, d, rank;
        slong *P;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        for (r = 0; r <= FLINT_MIN(m, n); r++)
        {
            TEMPLATE(T, mat_init) (A, m, n, ctx);
            TEMPLATE(T, mat_randrank) (A, state, r, ctx);

            if (n_randint(state, 2))
            {
                d = n_randint(state, 2 * m * n + 1);
                TEMPLATE(T, mat_randops) (A, d, state, ctx);
            }

            TEMPLATE(T, mat_init_set) (LU, A, ctx);
            P = flint_malloc(sizeof(slong) * m);

            rank = TEMPLATE(T, mat_lu_recursive) (P, LU, 0, ctx);

            if (r != rank)
            {
                printf("FAIL:\n");
                printf("wrong rank!\n");
                printf("A:");
                TEMPLATE(T, mat_print_pretty) (A, ctx);
                printf("LU:");
                TEMPLATE(T, mat_print_pretty) (LU, ctx);
                fflush(stdout);
                flint_abort();
            }

            check(P, LU, A, rank, ctx);

            TEMPLATE(T, mat_clear) (A, ctx);
            TEMPLATE(T, mat_clear) (LU, ctx);
            flint_free(P);
        }

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
