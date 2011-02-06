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
#include "fmpz_vec.h"
#include "fmpz_mat.h"



void _fmpz_mat_inv_2x2(fmpz ** b, fmpz_t den, fmpz ** const a)
{
    fmpz_t tmp;

    _fmpz_mat_det_cofactor_2x2(den, a);

    fmpz_neg(&b[0][1], &a[0][1]);
    fmpz_neg(&b[1][0], &a[1][0]);

    fmpz_init(tmp);
    fmpz_set(tmp, &a[0][0]);
    fmpz_set(&b[0][0], &a[1][1]);
    fmpz_set(&b[1][1], tmp);
    fmpz_clear(tmp);
}

void fmpz_mat_inv(fmpz_mat_t B, fmpz_t den, const fmpz_mat_t A)
{
    fmpz_mat_t I;
    long i, dim;

    dim = A->r;

    FMPZ_MAT_ASSERT(dim == A->c,
        "fmpz_mat_inv: matrix must be square");
    FMPZ_MAT_ASSERT(B->r == A->r,
        "fmpz_mat_inv: output must have same size as input");
    FMPZ_MAT_ASSERT(B->c == A->c,
        "fmpz_mat_inv: output must have same size as input");

    switch (dim)
    {
        case 0:
            fmpz_set_ui(den, 1UL);
            break;
        case 1:
            fmpz_set(den, A->entries);
            fmpz_set_ui(B->entries, 1UL);
            break;
        case 2:
            _fmpz_mat_inv_2x2(B->rows, den, A->rows);
            break;
        default:
            fmpz_mat_init(I, dim, dim);
            for (i = 0; i < dim; i++)
                fmpz_set_ui(&I->rows[i][i], 1UL);
            fmpz_mat_solve_mat(B, den, A, I);
            fmpz_mat_clear(I);
    }
}
