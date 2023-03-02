/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void fmpz_mat_mul_fmpz_vec(
    fmpz * c,
    const fmpz_mat_t A,
    const fmpz * b, slong blen)
{
    slong i;
    slong len = FLINT_MIN(A->c, blen);

    for (i = A->r - 1; i >= 0; i--)
        _fmpz_vec_dot(c + i, A->rows[i], b, len);
}

