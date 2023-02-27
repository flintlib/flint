/*
    Copyright (C) 2010,2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void
fmpz_mat_mul_classical(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
{
    slong ar, bc, br;
    slong i, j, k;

    ar = A->r;
    br = B->r;
    bc = B->c;

    if (br == 0)
    {
        fmpz_mat_zero(C);
        return;
    }

    for (i = 0; i < ar; i++)
    {
        for (j = 0; j < bc; j++)
        {
            fmpz_mul(fmpz_mat_entry(C, i, j),
                     fmpz_mat_entry(A, i, 0),
                     fmpz_mat_entry(B, 0, j));

            for (k = 1; k < br; k++)
            {
                fmpz_addmul(fmpz_mat_entry(C, i, j),
                            fmpz_mat_entry(A, i, k),
                            fmpz_mat_entry(B, k, j));
            }
        }
    }
}
