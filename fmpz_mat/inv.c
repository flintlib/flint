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

    Copyright (C) 2010-2011 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"

void
_fmpz_mat_inv_2x2(fmpz ** b, fmpz_t den, fmpz ** const a)
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

int
fmpz_mat_inv(fmpz_mat_t B, fmpz_t den, const fmpz_mat_t A)
{
    len_t dim = A->r;

    if (dim == 0)
    {
        fmpz_one(den);
        return 1;
    }
    else if (dim == 1)
    {
        fmpz_set(den, A->entries);
        fmpz_one(B->entries);
        return !fmpz_is_zero(den);
    }
    else if (dim == 2)
    {
        _fmpz_mat_inv_2x2(B->rows, den, A->rows);
        return !fmpz_is_zero(den);
    }
    else
    {
        fmpz_mat_t I;
        len_t i;
        int success;

        fmpz_mat_init(I, dim, dim);
        for (i = 0; i < dim; i++)
            fmpz_one(fmpz_mat_entry(I, i, i));
        success = fmpz_mat_solve_fflu(B, den, A, I);
        fmpz_mat_clear(I);
        return success;
    }
}
