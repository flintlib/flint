/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

int
fmpz_mat_col_equal(const fmpz_mat_t M, slong m, slong n)
{
   slong i;

   for (i = 0; i < M->r; i++)
   {
      if (!fmpz_equal(M->rows[i] + m, M->rows[i] + n))
         return 0;
   }

   return 1;
}
