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

    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include "fmpz_poly_mat.h"

void
fmpz_poly_mat_transpose(fmpz_poly_mat_t B, const fmpz_poly_mat_t A)
{
    slong i, j;

    if (B->r != A->c || B->c != A->r)
    {
        flint_printf("Exception (fmpz_poly_mat_transpose). Incompatible dimensions.\n");
        flint_abort();
    }

    if (A == B)  /* In-place, guaranteed to be square */
    {
        for (i = 0; i < A->r - 1; i++)
            for (j = i + 1; j < A->c; j++)
                fmpz_poly_swap(fmpz_poly_mat_entry(B, i, j), 
                               fmpz_poly_mat_entry(B, j, i));
    }
    else  /* Not aliased; general case */
    {
        for (i = 0; i < B->r; i++)
            for (j = 0; j < B->c; j++)
                fmpz_poly_set(fmpz_poly_mat_entry(B, i, j), 
                              fmpz_poly_mat_entry(A, j, i));
    }
}
