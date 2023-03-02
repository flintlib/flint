/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void
fmpz_mat_scalar_mul_2exp(fmpz_mat_t B, const fmpz_mat_t A, ulong exp)
{
    slong i, j;

    if (exp == 0)
    {
        fmpz_mat_set(B, A);
        return;
    }
    else if (exp == 1)
    {
        fmpz_mat_add(B, A, A);
        return;
    }

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            fmpz_mul_2exp(fmpz_mat_entry(B, i, j), fmpz_mat_entry(A, i, j),
                          exp);
}
