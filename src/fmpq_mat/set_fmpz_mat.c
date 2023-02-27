/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mat.h"

void fmpq_mat_set_fmpz_mat(fmpq_mat_t dest, const fmpz_mat_t src)
{
    slong i, j;

    for (i = 0; i < src->r; i++)
    {
        for (j = 0; j < src->c; j++)
        {
            fmpz_set(fmpq_mat_entry_num(dest, i, j), fmpz_mat_entry(src, i, j));
            fmpz_one(fmpq_mat_entry_den(dest, i, j));
        }
    }
}
