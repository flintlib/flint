/*
    Copyright (C) 2010-2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

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
    slong dim = A->r;

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
        slong i;
        int success;

        fmpz_mat_init(I, dim, dim);
        for (i = 0; i < dim; i++)
            fmpz_one(fmpz_mat_entry(I, i, i));
        success = fmpz_mat_solve_fflu(B, den, A, I);
        fmpz_mat_clear(I);
        return success;
    }
}
