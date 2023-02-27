/*
    Copyright (C) 2015 Dharak Kharod

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "fmpz.h"

void fmpz_mat_content(fmpz_t mat_gcd, const fmpz_mat_t A)
{
	slong rows = A->r;
	slong cols = A->c;
	slong i, j;

	fmpz_set_si(mat_gcd,0);

	for (i = 0; i < rows; i++ )
	{
		for (j = 0; j < cols; j++)
		{
			fmpz_gcd(mat_gcd,mat_gcd, fmpz_mat_entry(A, i, j));
			
			if (fmpz_is_one(mat_gcd))
			{
				return;
			}
		}
	}

	return;
}

