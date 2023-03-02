/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void fmpz_mat_scalar_smod(fmpz_mat_t B, const fmpz_mat_t A, const fmpz_t P)
{
   slong i;

   for (i = 0; i < A->r; i++)
      _fmpz_vec_scalar_smod_fmpz(B->rows[i], A->rows[i], A->c, P);
}
