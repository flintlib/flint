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
_fmpz_mat_mul(fmpz ** C, fmpz ** const A, long ar, long ac,
                         fmpz ** const B, long br, long bc)
{
    long i, j, k;

    for (i = 0; i < ar; i++)
    {
        for (j = 0; j < bc; j++)
        {
            fmpz_zero(&C[i][j]);
            for (k = 0; k < br; k++)
                fmpz_addmul(&C[i][j], &A[i][k], &B[k][j]);
        }
    }
}

void
fmpz_mat_mul(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
{
    long cr, cc;

    cr = A->r;
    cc = B->c;

    FMPZ_MAT_ASSERT((A->c == B->r && C->r == cr && C->c == cc),
        "fmpz_mat_mul: incompatible dimensions\n");

    if (C == A || C == B)
    {
        fmpz_mat_t t;
        fmpz_mat_init(t, cr, cc);
        _fmpz_mat_mul(t->rows, A->rows, A->r, A->c, B->rows, B->r, B->c);
        fmpz_mat_swap(C, t);
        fmpz_mat_clear(t);
    }
    else
    {
        _fmpz_mat_mul(C->rows, A->rows, A->r, A->c, B->rows, B->r, B->c);
    }
}
