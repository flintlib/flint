/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

slong
fmpz_mat_max_bits(const fmpz_mat_t mat)
{
    slong i;
    slong bits, row_bits, sign;

    sign = 1;
    bits = 0;

    if (mat->r == 0 || mat->c == 0)
        return 0;

    for (i = 0; i < mat->r; i++)
    {
        row_bits = _fmpz_vec_max_bits(mat->rows[i], mat->c);
        if (row_bits < 0)
        {
            row_bits = -row_bits;
            sign = -1;
        }
        bits = FLINT_MAX(bits, row_bits);
    }

    return bits * sign;
}
