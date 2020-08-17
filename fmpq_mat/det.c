/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mat.h"

void fmpq_mat_det(fmpq_t det, const fmpq_mat_t mat)
{
    slong n = mat->r;

    if (n == 0)
    {
        fmpq_set_si(det, WORD(1), WORD(1));
        return;
    }
    else if (n == 1)
    {
        fmpq_set(det, fmpq_mat_entry(mat, 0, 0));
    }
    else if (n == 2)
    {
        fmpq_t t;
        fmpq_init(t);

        fmpq_mul(t, fmpq_mat_entry(mat, 0, 0), fmpq_mat_entry(mat, 1, 1));
        fmpq_submul(t, fmpq_mat_entry(mat, 0, 1), fmpq_mat_entry(mat, 1, 0));

        fmpq_set(det, t);
        fmpq_clear(t);
    }
    else
    {
        fmpz_mat_t num;
        fmpz * den;
        slong i;

        fmpz_mat_init(num, mat->r, mat->c);
        den = _fmpz_vec_init(mat->r);

        fmpq_mat_get_fmpz_mat_rowwise(num, den, mat);
        fmpz_mat_det(&det->num, num);

        fmpz_one(&det->den);
        for (i = 0; i < mat->r; i++)
            fmpz_mul(&det->den, &det->den, den + i);

        fmpq_canonicalise(det);

        fmpz_mat_clear(num);
        _fmpz_vec_clear(den, mat->r);
    }
}
