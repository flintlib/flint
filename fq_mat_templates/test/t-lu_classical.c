/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"

#include <stdio.h>
#include <stdlib.h>
#include "ulong_extras.h"

void perm(TEMPLATE(T, mat_t) A, slong * P)
{
    slong i;
    TEMPLATE(T, struct) ** tmp;

    if (A->c == 0 || A->r == 0)
        return;

    tmp = flint_malloc(sizeof(TEMPLATE(T, struct)*) * A->r);

    for (i = 0; i < A->r; i++) tmp[P[i]] = A->rows[i];
    for (i = 0; i < A->r; i++) A->rows[i] = tmp[i];

    flint_free(tmp);
}

void check(slong * P, TEMPLATE(T, mat_t) LU, const TEMPLATE(T, mat_t) A, slong rank,
           const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, mat_t) B, L, U;
    slong m, n, i, j;

    m = A->r;
    n = A->c;

    TEMPLATE(T, mat_init)(B, m, n, ctx);
    TEMPLATE(T, mat_init)(L, m, m, ctx);
    TEMPLATE(T, mat_init)(U, m, n, ctx);

    rank = FLINT_ABS(rank);

    for (i = rank; i < FLINT_MIN(m, n); i++)
    {
        for (j = i; j < n; j++)
        {
            if (!TEMPLATE(T, is_zero)(TEMPLATE(T, mat_entry)(LU, i, j), ctx))
            {
                printf("FAIL: wrong shape!\n");
                abort();
            }
        }
    }

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < FLINT_MIN(i, n); j++)
            TEMPLATE(T, mat_entry_set)(L, i, j, TEMPLATE(T, mat_entry)(LU, i, j));
        if (i < rank)
            TEMPLATE(T, one)(TEMPLATE(T, mat_entry)(L, i, i), ctx);
        for (j = i; j < n; j++)
            TEMPLATE(T, mat_entry_set)(U, i, j, TEMPLATE(T, mat_entry)(LU, i, j));
    }

    TEMPLATE(T, mat_mul)(B, L, U, ctx);
    perm(B, P);

    if (!TEMPLATE(T, mat_equal)(A, B))
    {
        printf("FAIL\n");
        printf("A:\n");
        TEMPLATE(T, mat_print_pretty)(A, ctx);
        printf("LU:\n");
        TEMPLATE(T, mat_print_pretty)(LU, ctx);
        printf("B:\n");
        TEMPLATE(T, mat_print_pretty)(B, ctx);
        abort();
    }

    TEMPLATE(T, mat_clear)(B);
    TEMPLATE(T, mat_clear)(L);
    TEMPLATE(T, mat_clear)(U);
}



int
main(void)
{
    slong i;

    flint_rand_t state;
    flint_randinit(state);

    printf("lu_classical....");
    fflush(stdout);

    for (i = 0; i < 100; i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, mat_t) A, LU;
        
        slong m, n, r, d, rank;
        slong * P;

        TEMPLATE(T, ctx_randtest)(ctx, state);
        
        m = n_randint(state, 20);
        n = n_randint(state, 20);

        for (r = 0; r <= FLINT_MIN(m, n); r++)
        {
            TEMPLATE(T, mat_init)(A, m, n, ctx);
            TEMPLATE(T, mat_randrank)(A, state, r, ctx);

            if (n_randint(state, 2))
            {
                d = n_randint(state, 2*m*n + 1);
                TEMPLATE(T, mat_randops)(A, d, state, ctx);
            }

            TEMPLATE(T, mat_init_set)(LU, A, ctx);
            P = flint_malloc(sizeof(slong) * m);

            rank = TEMPLATE(T, mat_lu_classical)(P, LU, 0, ctx);

            if (r != rank)
            {
                printf("FAIL:\n");
                printf("wrong rank!\n");
                printf("A:");
                TEMPLATE(T, mat_print_pretty)(A, ctx);
                printf("LU:");
                TEMPLATE(T, mat_print_pretty)(LU, ctx);
                abort();
            }

            check(P, LU, A, rank, ctx);

            TEMPLATE(T, mat_clear)(A);
            TEMPLATE(T, mat_clear)(LU);
            flint_free(P);

        }

        TEMPLATE(T, ctx_clear)(ctx);
    }

    flint_randclear(state);
    printf("PASS\n");
    return 0;
}


#endif
