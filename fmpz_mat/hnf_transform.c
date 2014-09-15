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

    Copyright (C) 2014 Alex J. Best

******************************************************************************/

#include "fmpz_mat.h"

void
fmpz_mat_hnf_transform(fmpz_mat_t H, fmpz_mat_t U, const fmpz_mat_t A)
{
    slong i, j, m, n;
    fmpz_mat_t A2, H2;

    m = fmpz_mat_nrows(A);
    n = fmpz_mat_ncols(A);

    fmpz_mat_init(A2, m, n + m);
    fmpz_mat_init(H2, m, n + m);

    /* append an identity matrix to the right of A */
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
            fmpz_set(fmpz_mat_entry(A2, i, j), fmpz_mat_entry(A, i, j));
        fmpz_one(fmpz_mat_entry(A2, i, n + i));
    }

    fmpz_mat_hnf(H2, A2);

    /* recover H and U */
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
            fmpz_set(fmpz_mat_entry(H, i, j), fmpz_mat_entry(H2, i, j));
        for (j = n; j < n + m; j++)
            fmpz_set(fmpz_mat_entry(U, i, j - n), fmpz_mat_entry(H2, i, j));
    }

    fmpz_mat_clear(A2);
    fmpz_mat_clear(H2);
}
