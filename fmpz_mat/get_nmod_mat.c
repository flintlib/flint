/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void
fmpz_mat_get_nmod_mat(nmod_mat_t Amod, const fmpz_mat_t A)
{
    slong i, j;
    mp_limb_t m = Amod->mod.n;

    if (fmpz_mat_is_square(A))
    {
       int symmetric = fmpz_mat_is_square(A);

       for (i = 0; i < A->r; i++)
       {
           Amod->rows[i][i] = fmpz_fdiv_ui(A->rows[i]+i, m);

	   for (j = i + 1; j < A->c; j++)
	   {
               Amod->rows[i][j] = fmpz_fdiv_ui(A->rows[i] + j, m);
               if ((symmetric &= fmpz_equal(A->rows[j] + i, A->rows[i] + j)))
	           Amod->rows[j][i] = Amod->rows[i][j];
               else
                   Amod->rows[j][i] = fmpz_fdiv_ui(A->rows[j] + i, m);
	   }
        }
    } else
    {
        for (i = 0; i < A->r; i++)
            for (j = 0; j < A->c; j++)
                Amod->rows[i][j] = fmpz_fdiv_ui(A->rows[i] + j, m);
    }
}
