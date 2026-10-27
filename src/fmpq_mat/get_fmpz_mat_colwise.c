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
fmpq_mat_get_fmpz_mat_colwise(fmpz_mat_t num, fmpz * den, const fmpq_mat_t mat)
{
    fmpz_t lcm;
    slong i, j;

    if (fmpq_mat_is_empty(mat))
        return;

    fmpz_init(lcm);

    for (j = 0; j < mat->c; j++)
    {
        /* Compute common denominator of column (the denominators are
           2 stride fmpz apart) */
        _fmpz_vec_lcm_strided(lcm, fmpq_mat_entry_den(mat, 0, j), 2 * mat->stride, mat->r, NULL);

        if (den != NULL)
            fmpz_set(den + j, lcm);

        /* Rescale numerators in column */
        if (fmpz_is_one(lcm))
        {
            for (i = 0; i < mat->r; i++)
                fmpz_set(fmpz_mat_entry(num, i, j),
                         fmpq_mat_entry_num(mat, i, j));
        }
        else
        {
            _fmpq_vec_get_fmpz_vec_scaled(fmpz_mat_entry(num, 0, j), num->stride,
                fmpq_mat_entry(mat, 0, j), mat->stride, mat->r, lcm);
        }
    }

    fmpz_clear(lcm);
}
