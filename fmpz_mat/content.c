/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

  Copyright (C) 2015 Dharak Kharod

******************************************************************************/

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

