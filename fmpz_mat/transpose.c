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

******************************************************************************/

#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"

void
fmpz_mat_transpose(fmpz_mat_t B, const fmpz_mat_t A)
{
    fmpz tmp;
    len_t i, j;

    if (B->r != A->c || B->c != A->r)
    {
        printf("Exception (fmpz_mat_transpose). Incompatible dimensions.\n");
        abort();
    }

    if (A == B)  /* In-place, guaranteed to be square */
    {
        for (i = 0; i < A->r - 1; i++)
            for (j = i + 1; j < A->c; j++)
            {
                tmp = A->rows[i][j];
                A->rows[i][j] = A->rows[j][i];
                A->rows[j][i] = tmp;
            }
    }
    else  /* Not aliased; general case */
    {
        for (i = 0; i < B->r; i++)
            for (j = 0; j < B->c; j++)
                fmpz_set(&B->rows[i][j], &A->rows[j][i]);
    }
}
