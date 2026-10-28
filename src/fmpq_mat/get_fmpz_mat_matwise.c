/*
    Copyright (C) 2011, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_vec.h"
#include "fmpq_vec.h"
#include "fmpq_mat.h"

void
fmpq_mat_get_fmpz_mat_matwise(fmpz_mat_t num, fmpz_t den, const fmpq_mat_t mat)
{
    fmpz_t lcm;
    slong i, j;

    if (fmpq_mat_is_empty(mat))
    {
        fmpz_one(den);
        return;
    }

    fmpz_init(lcm);

    /* Compute common denominator of matrix, row by row */
    _fmpz_vec_lcm_strided(lcm, fmpq_mat_entry_den(mat, 0, 0), 2, mat->c, NULL);
    for (i = 1; i < mat->r; i++)
        _fmpz_vec_lcm_strided(lcm, fmpq_mat_entry_den(mat, i, 0), 2, mat->c, lcm);

    fmpz_set(den, lcm);

    for (i = 0; i < mat->r; i++)
    {
        /* Rescale numerators */
        if (fmpz_is_one(lcm))
        {
            for (j = 0; j < mat->c; j++)
                fmpz_set(fmpz_mat_entry(num, i, j),
                         fmpq_mat_entry_num(mat, i, j));
        }
        else
        {
            _fmpq_vec_get_fmpz_vec_scaled(fmpz_mat_entry(num, i, 0), 1,
                fmpq_mat_entry(mat, i, 0), 1, mat->c, lcm);
        }
    }

    fmpz_clear(lcm);
}
