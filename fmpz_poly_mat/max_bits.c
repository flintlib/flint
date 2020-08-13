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
fmpz_poly_mat_max_bits(const fmpz_poly_mat_t A)
{
    slong i, j, bits, max;
    int sign;

    max = 0; 
    sign = 0;

    for (i = 0; i < A->r; i++)
    {
        for (j = 0; j < A->c; j++)
        {
            bits = fmpz_poly_max_bits(fmpz_poly_mat_entry(A, i, j));

            if (bits < 0)
            {
                sign = 1;
                max = FLINT_MAX(max, -bits);
            }
            else
                max = FLINT_MAX(max, bits);
        }
    }

    return sign ? -max : max;
}
