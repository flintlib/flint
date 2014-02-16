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

    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include "fmpz_mat.h"

void fmpz_mat_gram(fmpz_mat_t B, const fmpz_mat_t A)
{
	slong i, j, k;
	
	if(B->r != A->r || B->c != A->r) {
		flint_printf("Exception (fmpz_mat_gram). Incompatible dimensions.\n");
		abort();
	}
	
	if(B == A) {
		fmpz_mat_t t;
		fmpz_mat_init(t, B->r, B->c);
		fmpz_mat_gram(t, A);
		fmpz_mat_swap(B, t);
		fmpz_mat_clear(t);
		return;
	}
	
	if(A->c == 0) {
		fmpz_mat_zero(B);
		return;
	}
	
	for(i = 0; i < B->r; i++) {
		for(j = 0; j < B->c; j++) {
			fmpz_mul(fmpz_mat_entry(B, i, j),
					 fmpz_mat_entry(A, i, 0),
					 fmpz_mat_entry(A, j, 0));
					 
			for (k = 1; k < A->c; k++) {
                fmpz_addmul(fmpz_mat_entry(B, i, j),
                            fmpz_mat_entry(A, i, k),
                            fmpz_mat_entry(A, j, k));
            }
		}
	}
}
