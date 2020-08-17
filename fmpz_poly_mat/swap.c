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

void
fmpz_poly_mat_swap(fmpz_poly_mat_t A, fmpz_poly_mat_t B)
{
    if (A != B)
    {
        fmpz_poly_mat_struct tmp;

        tmp = *A;
        *A = *B;
        *B = tmp;
    }
}
