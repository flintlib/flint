/*
    Copyright (C) 2015 Anubhav Srivastava
    Copyright (C) 2015 Elena Sergeicheva

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_mat.h"

void
fmpz_poly_mat_concat_horizontal(fmpz_poly_mat_t res, const fmpz_poly_mat_t mat1, const fmpz_poly_mat_t mat2)
{
    slong i, j;
    slong r1 = mat1->r;
    slong c1 = mat1->c;
    slong r2 = mat2->r;
    slong c2 = mat2->c;
    
    for (i = 0; i < r1; i++)
    {
        for (j = 0; j < c1; j++)
        {
            fmpz_poly_set(fmpz_poly_mat_entry(res, i, j), fmpz_poly_mat_entry(mat1, i, j));
        }
    }

    for (i = 0; i < r2; i++)
    {
        for (j = 0; j < c2; j++)
        {
            fmpz_poly_set(fmpz_poly_mat_entry(res, i, j + c1), fmpz_poly_mat_entry(mat2, i, j));
        }
    }
}
