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
#include "flint.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"

int
fmpz_poly_mat_pivot(long * perm, fmpz_poly_mat_t A, long r, long c)
{
    long t, j;
    fmpz_poly_struct * u;

    if (!fmpz_poly_is_zero(fmpz_poly_mat_entry(A, r, c)))
        return 1;

    for (j = r + 1; j < A->r; j++)
    {
        if (!fmpz_poly_is_zero(fmpz_poly_mat_entry(A, j, c)))
        {
            if (perm)
            {
                t = perm[j];
                perm[j] = perm[r];
                perm[r] = t;
            }

            u = A->rows[j];
            A->rows[j] = A->rows[r];
            A->rows[r] = u; 
            return -1;
        }
    }
    return 0;
}
