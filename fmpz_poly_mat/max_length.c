/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_mat-impl.h"


slong
fmpz_poly_mat_max_length(const fmpz_poly_mat_t A)
{
    slong i, j, len, max;

    max = 0;

    for (i = 0; i < A->r; i++)
    {
        for (j = 0; j < A->c; j++)
        {
            len = fmpz_poly_length(fmpz_poly_mat_entry(A, i, j));
            max = FLINT_MAX(len, max);
        }
    }

    return max;
}
