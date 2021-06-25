/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void fmpz_mat_fmpz_vec_mul_ptr(
    fmpz * const * c,
    const fmpz * const * a, slong alen,
    const fmpz_mat_t B)
{
    slong i, j;
    slong len = FLINT_MIN(B->r, alen);

    for (i = B->c - 1; i >= 0; i--)
    {
        fmpz_zero(c[i]);
        for (j = 0; j < len; j++)
            fmpz_addmul(c[i], a[j], fmpz_mat_entry(B, j, i));
    }
}

