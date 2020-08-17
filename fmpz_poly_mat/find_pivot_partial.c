/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "flint.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"

slong
fmpz_poly_mat_find_pivot_partial(const fmpz_poly_mat_t mat,
                                    slong start_row, slong end_row, slong c)
{
    slong best_row, best_length, best_bits, i;

    best_row = start_row;
    best_length = fmpz_poly_length(fmpz_poly_mat_entry(mat, start_row, c));

    best_bits = fmpz_poly_max_bits(fmpz_poly_mat_entry(mat, start_row, c));
    best_bits = FLINT_ABS(best_bits);

    for (i = start_row + 1; i < end_row; i++)
    {
        slong b, l;

        l = fmpz_poly_length(fmpz_poly_mat_entry(mat, i, c));

        if (l != 0 && (best_length == 0 || l <= best_length))
        {
            b = fmpz_poly_max_bits(fmpz_poly_mat_entry(mat, i, c));
            b = FLINT_ABS(b);

            if (best_length == 0 || l < best_length || b < best_bits)
            {
                best_row = i;
                best_length = l;
                best_bits = b;
            }
        }
    }

    if (best_length == 0)
        return -1;

    return best_row;
}
