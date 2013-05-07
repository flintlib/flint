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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"
#include "perm.h"

long
nmod_mat_rref(nmod_mat_t A)
{
    long i, j, k, n, rank;
    long * pivots;
    long * nonpivots;
    long * P;

    nmod_mat_t U, V;

    n = A->c;

    P = _perm_init(nmod_mat_nrows(A));
    rank = nmod_mat_lu(P, A, 0);
    _perm_clear(P);

    if (rank == 0)
        return rank;

    /* Clear L */
    for (i = 0; i < A->r; i++)
        for (j = 0; j < FLINT_MIN(i, rank); j++)
            nmod_mat_entry(A, i, j) = 0UL;

    /* We now reorder U to proper upper triangular form U | V
       with U full-rank triangular, set V = U^(-1) V, and then
       put the column back in the original order.

       An improvement for some matrices would be to compress V by
       discarding columns containing nothing but zeros. */

    nmod_mat_init(U, rank, rank, A->mod.n);
    nmod_mat_init(V, rank, n - rank, A->mod.n);

    pivots = flint_malloc(sizeof(long) * rank);
    nonpivots = flint_malloc(sizeof(long) * (n - rank));

    for (i = j = k = 0; i < rank; i++)
    {
        while (nmod_mat_entry(A, i, j) == 0UL)
        {
            nonpivots[k] = j;
            k++;
            j++;
        }
        pivots[i] = j;
        j++;
    }
    while (k < n - rank)
    {
        nonpivots[k] = j;
        k++;
        j++;
    }

    for (i = 0; i < rank; i++)
    {
        for (j = 0; j <= i; j++)
            nmod_mat_entry(U, j, i) = nmod_mat_entry(A, j, pivots[i]);
    }

    for (i = 0; i < n - rank; i++)
    {
        for (j = 0; j < rank; j++)
            nmod_mat_entry(V, j, i) = nmod_mat_entry(A, j, nonpivots[i]);
    }

    nmod_mat_solve_triu(V, U, V, 0);

    /* Clear pivot columns */
    for (i = 0; i < rank; i++)
    {
        for (j = 0; j <= i; j++)
            nmod_mat_entry(A, j, pivots[i]) = (i == j);
    }

    /* Write back the actual content */
    for (i = 0; i < n - rank; i++)
    {
        for (j = 0; j < rank; j++)
            nmod_mat_entry(A, j, nonpivots[i]) = nmod_mat_entry(V, j, i);
    }

    nmod_mat_clear(U);
    nmod_mat_clear(V);

    flint_free(pivots);
    flint_free(nonpivots);

    return rank;
}
