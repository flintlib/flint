/*
    Copyright (C) 2015 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

slong
fmpz_mat_howell_form_mod(fmpz_mat_t A, const fmpz_t mod)
{
    slong i, j, n;
    slong k;

    if (fmpz_mat_is_empty(A))
        return 0;

    n = A->r;
    k = n;
    fmpz_mat_strong_echelon_form_mod(A, mod);

    for (i = 0; i < n; i++)
    {
        if (fmpz_mat_is_zero_row(A, i))
        {
            k--;
            for (j = i + 1; j < n; j++)
            {
                if (!fmpz_mat_is_zero_row(A, j))
                {
                    fmpz_mat_swap_rows(A, NULL, i, j);
                    j = n;
                    k++;
                }
            }
        }
    }
    return k;
}

