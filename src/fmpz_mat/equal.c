/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"

int fmpz_mat_equal(const fmpz_mat_t mat1, const fmpz_mat_t mat2)
{
    slong j;

    if (mat1->r != mat2->r || mat1->c != mat2->c)
    {
        return 0;
    }

    if (mat1->r == 0 || mat1->c == 0)
        return 1;

    for (j = 0; j < mat1->r; j++)
    {
        if (!_fmpz_vec_equal(fmpz_mat_row(mat1, j), fmpz_mat_row(mat2, j), mat1->c))
        {
            return 0;
        }
    }

    return 1;
}

int fmpz_mat_equal_col(fmpz_mat_t M, slong m, slong n)
{
   slong i;

   for (i = 0; i < M->r; i++)
   {
      if (!fmpz_equal(fmpz_mat_entry(M, i, m), fmpz_mat_entry(M, i, n)))
         return 0;
   }

   return 1;
}

int fmpz_mat_equal_row(fmpz_mat_t M, slong m, slong n)
{
   slong i;

   for (i = 0; i < M->c; i++)
   {
      if (!fmpz_equal(fmpz_mat_entry(M, m, i), fmpz_mat_entry(M, n, i)))
         return 0;
   }

   return 1;
}
