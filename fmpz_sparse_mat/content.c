/*
    Copyright (C) 2015 Dharak Kharod

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_sparse_mat.h"
#include "fmpz.h"

void fmpz_sparse_mat_content(fmpz_t mat_gcd, const fmpz_sparse_mat_t M)
{
	slong i, j;
	fmpz_set_si(mat_gcd,0);
	for (i = 0; i < M->r; i++ )
	{
		for (j = 0; j < M->rows[i].nnz; j++)
		{
			fmpz_gcd(mat_gcd, mat_gcd, M->rows[i].entries[j].val);
			if (fmpz_is_one(mat_gcd)) return;
		}
	}
}

