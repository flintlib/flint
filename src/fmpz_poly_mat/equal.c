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

int
fmpz_poly_mat_equal(const fmpz_poly_mat_t A, const fmpz_poly_mat_t B)
{
    slong i, j;

    if (A->r != B->r || A->c != B->c)
        return 0;

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            if (!fmpz_poly_equal(fmpz_poly_mat_entry(A, i, j),
                                 fmpz_poly_mat_entry(B, i, j)))
                return 0;
    return 1;
}
