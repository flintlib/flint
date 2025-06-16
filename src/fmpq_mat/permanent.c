/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpq.h"
#include "fmpq_mat.h"

int
fmpq_mat_permanent(fmpq_t res, const fmpq_mat_t mat)
{
    slong n = mat->r;

    if (n == 0)
    {
        fmpq_one(res);
        return 1;
    }
    else if (n == 1)
    {
        fmpq_set(res, fmpq_mat_entry(mat, 0, 0));
        return 1;
    }
    else if (n == 2)
    {
        fmpq_mul(res, fmpq_mat_entry(mat, 0, 0), fmpq_mat_entry(mat, 1, 1));
        fmpq_addmul(res, fmpq_mat_entry(mat, 0, 1), fmpq_mat_entry(mat, 1, 0));
        return 1;
    }
    else
    {
        fmpz_mat_t num;
        fmpz * den;
        slong i;
        int success;

        fmpz_mat_init(num, mat->r, mat->c);
        den = _fmpz_vec_init(mat->r);

        fmpq_mat_get_fmpz_mat_rowwise(num, den, mat);
        success = fmpz_mat_permanent(fmpq_numref(res), num);

        fmpz_mul(fmpq_denref(res), den, den + 1);
        for (i = 2; i < mat->r; i++)
            fmpz_mul(fmpq_denref(res), fmpq_denref(res), den + i);

        fmpq_canonicalise(res);

        fmpz_mat_clear(num);
        _fmpz_vec_clear(den, mat->r);

        return success;
    }
}

