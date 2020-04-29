/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_sparse_mat.h"

void fmpz_sparse_mat_gram(fmpz_mat_t B, const fmpz_sparse_mat_t M)
{
	slong i, j;

	if (M->r == 0) fmpz_mat_zero(B);	
	else
		for (i = 0; i < M->r; i++) 
			for (j = 0; j < M->r; j++) 
				fmpz_sparse_vec_dot(&B->rows[i][j], &M->rows[i], &M->rows[j]);
}
