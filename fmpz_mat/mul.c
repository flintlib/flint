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
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"

void
_fmpz_mat_mul(fmpz * C,
              const fmpz * A, long ar, long ac,
              const fmpz * B, long br, long bc)
{
    long i, j, k;
    fmpz * a, b, c;

    for (i = 0; i < ar; i++)
    {
        for (j = 0; j < bc; j++)
        {
            a = A + ac*i;
            c = C + bc*i + j;
            fmpz_zero(c);

            for (k = 0, b = B + j; k < br; k++, a++)
            {
                /* C[i,j] += A[i,k] * B[k,j] */
                b = B + j + bc*k;
                fmpz_addmul(c, a, b);
            }
        }
    }
}

void
fmpz_mat_mul(fmpz_mat_t C, fmpz_mat_t A, fmpz_mat_t B)
{
    if (A->c != B->r || C->r != A->r || C->c != B->c)
    {
        printf("fmpz_mat_mul: incompatible dimensions\n");
        abort();
    }

    if (C == A || C == B)
    {
        printf("fmpz_mat_mul: aliasing not implemented\n");
        abort();
    }

    _fmpz_mat_mul(C->entries, A->entries, A->r, A->c, B->entries, B->r, B->c);
}
