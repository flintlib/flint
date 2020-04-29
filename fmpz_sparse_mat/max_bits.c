/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_sparse_mat.h"
#include "hashmap.h"
#include "longlong.h"

slong fmpz_sparse_mat_max_bits(const fmpz_sparse_mat_t M)
{
    slong i;
    slong bits, row_bits, sign;

    sign = 1;
    bits = 0;

    if (M->r == 0 || M->c == 0)
        return 0;

    for (i = 0; i < M->r; i++)
    {
        row_bits = fmpz_sparse_vec_max_bits(&M->rows[i]);
        if (row_bits < 0)
        {
            row_bits = -row_bits;
            sign = -1;
        }
        bits = FLINT_MAX(bits, row_bits);
    }
    return bits * sign;
}