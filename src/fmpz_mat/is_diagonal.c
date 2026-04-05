/*
    Copyright (C) 2023 Jean Kieffer
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"

int fmpz_mat_is_diagonal(const fmpz_mat_t A)
{
    slong r = fmpz_mat_nrows(A);
    slong c = fmpz_mat_ncols(A);
    slong j, k;

    for (j = 0; j < r; j++)
    {
        for (k = 0; k < c; k++)
        {
            if (j != k && !fmpz_is_zero(fmpz_mat_entry(A, j, k)))
            {
                return 0;
            }
        }
    }
    return 1;
}
