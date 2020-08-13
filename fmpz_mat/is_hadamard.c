/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

int
fmpz_mat_is_hadamard(const fmpz_mat_t A)
{
    slong i, j, n;
    fmpz_mat_t B, C;
    int result;

    n = fmpz_mat_nrows(A);

    if (n != fmpz_mat_ncols(A))
        return 0;

    if (n > 2 && n % 4 != 0)
        return 0;

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            if (!fmpz_is_pm1(fmpz_mat_entry(A, i, j)))
                return 0;

    fmpz_mat_init(B, n, n);
    fmpz_mat_init(C, n, n);

    fmpz_mat_transpose(B, A);
    fmpz_mat_mul(C, A, B);

    result = 1;

    for (i = 0; i < n && result; i++)
        for (j = 0; j < n && result; j++)
            result = (*fmpz_mat_entry(C, i, j) == n * (i == j));

    fmpz_mat_clear(B);
    fmpz_mat_clear(C);

    return result;
}

