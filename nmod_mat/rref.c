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
#include <mpir.h>
#include "flint.h"
#include "nmod_mat.h"

long
nmod_mat_rref(long * P, nmod_mat_t A)
{
    long i, j, col, rank;
    nmod_mat_t U, V;

    rank = nmod_mat_lu(P, A, 0);

    if (rank == 0)
        return rank;

    /* Clear L */
    for (i = 0; i < A->r; i++)
        for (j = 0; j < FLINT_MIN(i, rank); j++)
            nmod_mat_entry(A, i, j) = 0UL;

    nmod_mat_init(V, rank, rank, A->mod.n);
    nmod_mat_window_init(U, A, 0, 0, rank, A->c);

    /* V = compact(U) */
    col = 0;
    for (i = 0; i < rank; i++)
    {
        while (nmod_mat_entry(A, i, col) == 0UL)
            col++;
        for (j = 0; j <= i; j++)
            nmod_mat_entry(V, j, i) = nmod_mat_entry(A, j, col);
    }

    /* R = V^{-1} U */
    nmod_mat_solve_triu(U, V, U, 0);

    nmod_mat_clear(V);
    nmod_mat_window_clear(U);

    return rank;
}
