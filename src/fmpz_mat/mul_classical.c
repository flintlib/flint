/*
    Copyright (C) 2010, 2011, 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"

void
fmpz_mat_mul_classical(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
{
    slong ar, br, bc, i, j;

    ar = fmpz_mat_nrows(A);
    br = fmpz_mat_nrows(B);
    bc = fmpz_mat_ncols(B);

    if (ar == 0 || br == 0 || bc == 0)
    {
        fmpz_mat_zero(C);
        return;
    }

    if (br == 1)
    {
        for (i = 0; i < ar; i++)
        {
            for (j = 0; j < bc; j++)
            {
                fmpz_mul(fmpz_mat_entry(C, i, j),
                                 fmpz_mat_entry(A, i, 0),
                                 fmpz_mat_entry(B, 0, j));
            }
        }
    }
    else if (br == 2)
    {
        for (i = 0; i < ar; i++)
        {
            for (j = 0; j < bc; j++)
            {
                fmpz_fmma(fmpz_mat_entry(C, i, j),
                                 fmpz_mat_entry(A, i, 0),
                                 fmpz_mat_entry(B, 0, j),
                                 fmpz_mat_entry(A, i, 1),
                                 fmpz_mat_entry(B, 1, j));
            }
        }
    }
    else
    {
        fmpz * tmp;
        TMP_INIT;

        TMP_START;
        tmp = TMP_ALLOC(sizeof(fmpz) * br * bc);

        for (i = 0; i < br; i++)
            for (j = 0; j < bc; j++)
                tmp[j * br + i] = *fmpz_mat_entry(B, i, j);

        for (i = 0; i < ar; i++)
        {
            for (j = 0; j < bc; j++)
            {
                _fmpz_vec_dot_general(fmpz_mat_entry(C, i, j),
                    NULL, 0, fmpz_mat_entry(A, i, 0), tmp + j * br, 0, br);
            }
        }

        TMP_END;
    }
}
