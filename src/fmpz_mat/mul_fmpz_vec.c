/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"

void fmpz_mat_mul_fmpz_vec(
    fmpz * c,
    const fmpz_mat_t A,
    const fmpz * b, slong blen)
{
    slong i;
    slong len = FLINT_MIN(A->c, blen);

    for (i = A->r - 1; i >= 0; i--)
        _fmpz_vec_dot(c + i, fmpz_mat_row(A, i), b, len);
}

void fmpz_mat_mul_fmpz_vec_ptr(
    fmpz * const * c,
    const fmpz_mat_t A,
    const fmpz * const * b, slong blen)
{
    slong i, j;
    slong len = FLINT_MIN(A->c, blen);

    for (i = A->r - 1; i >= 0; i--)
    {
        fmpz * ci = c[i];
        const fmpz * Ai = fmpz_mat_row(A, i);

        fmpz_zero(ci);
        for (j = 0; j < len; j++)
            fmpz_addmul(ci, Ai + j, b[j]);
    }
}
