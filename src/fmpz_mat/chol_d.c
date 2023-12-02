/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"

#ifdef __GNUC__
# define sqrt __builtin_sqrt
#else
# include <math.h>
#endif

void
fmpz_mat_chol_d(d_mat_t R, const fmpz_mat_t A)
{
    slong i, k, j, r = A->r;

    if (A->r != A->c || R->r != A->r || R->c != A->c)
    {
        flint_throw(FLINT_ERROR, "(fmpz_mat_chol_d): Incompatible dimensions.\n");
    }

    if (A->r == 0)
    {
        return;
    }

    for (i = 0; i < r; i++)
    {
        for (j = 0; j < i + 1; j++)
        {
            double s = 0;
            for (k = 0; k < j; k++)
            {
                s += d_mat_entry(R, i, k) * d_mat_entry(R, j, k);
            }
            if (i == j)
                d_mat_entry(R, i, j) =
                    sqrt(fmpz_get_d(fmpz_mat_entry(A, i, i)) - s);
            else
                d_mat_entry(R, i, j) =
                    (fmpz_get_d(fmpz_mat_entry(A, i, j)) - s) / d_mat_entry(R,
                                                                            j,
                                                                            j);
        }
    }
}
