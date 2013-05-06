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

    Copyright (C) 2010,2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "ulong_extras.h"

void perm(nmod_mat_t A, long * P)
{
    long i;
    mp_ptr * tmp;

    if (A->c == 0 || A->r == 0)
        return;

    tmp = flint_malloc(sizeof(mp_ptr) * A->r);

    for (i = 0; i < A->r; i++) tmp[P[i]] = A->rows[i];
    for (i = 0; i < A->r; i++) A->rows[i] = tmp[i];

    flint_free(tmp);
}

void check(long * P, nmod_mat_t LU, const nmod_mat_t A, long rank)
{
    nmod_mat_t B, L, U;
    long m, n, i, j;

    m = A->r;
    n = A->c;

    nmod_mat_init(B, m, n, A->mod.n);
    nmod_mat_init(L, m, m, A->mod.n);
    nmod_mat_init(U, m, n, A->mod.n);

    rank = FLINT_ABS(rank);

    for (i = rank; i < FLINT_MIN(m, n); i++)
    {
        for (j = i; j < n; j++)
        {
            if (nmod_mat_entry(LU, i, j) != 0)
            {
                printf("FAIL: wrong shape!\n");
                abort();
            }
        }
    }

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < FLINT_MIN(i, n); j++)
            nmod_mat_entry(L, i, j) = nmod_mat_entry(LU, i, j);
        if (i < rank)
            nmod_mat_entry(L, i, i) = 1UL;
        for (j = i; j < n; j++)
            nmod_mat_entry(U, i, j) = nmod_mat_entry(LU, i, j);
    }

    nmod_mat_mul(B, L, U);
    perm(B, P);

    if (!nmod_mat_equal(A, B))
    {
        printf("FAIL\n");
        printf("A:\n");
        nmod_mat_print_pretty(A);
        printf("LU:\n");
        nmod_mat_print_pretty(LU);
        printf("B:\n");
        nmod_mat_print_pretty(B);
        abort();
    }

    nmod_mat_clear(B);
    nmod_mat_clear(L);
    nmod_mat_clear(U);
}



int
main(void)
{
    long i;

    flint_rand_t state;
    flint_randinit(state);

    printf("lu_classical....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, LU;
        mp_limb_t mod;
        long m, n, r, d, rank;
        long * P;

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        mod = n_randtest_prime(state, 0);

        for (r = 0; r <= FLINT_MIN(m, n); r++)
        {
            nmod_mat_init(A, m, n, mod);
            nmod_mat_randrank(A, state, r);

            if (n_randint(state, 2))
            {
                d = n_randint(state, 2*m*n + 1);
                nmod_mat_randops(A, d, state);
            }

            nmod_mat_init_set(LU, A);
            P = flint_malloc(sizeof(long) * m);

            rank = nmod_mat_lu_classical(P, LU, 0);

            if (r != rank)
            {
                printf("FAIL:\n");
                printf("wrong rank!\n");
                printf("A:");
                nmod_mat_print_pretty(A);
                printf("LU:");
                nmod_mat_print_pretty(LU);
                abort();
            }

            check(P, LU, A, rank);

            nmod_mat_clear(A);
            nmod_mat_clear(LU);
            flint_free(P);
        }
    }

    flint_randclear(state);
    printf("PASS\n");
    return 0;
}
