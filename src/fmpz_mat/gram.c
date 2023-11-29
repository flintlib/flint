/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"

void fmpz_mat_gram(fmpz_mat_t B, const fmpz_mat_t A)
{
	slong i, j, k;

	if (B->r != A->r || B->c != A->r)
    {
        flint_throw(FLINT_ERROR, "(fmpz_mat_gram): Incompatible dimensions.\n");
	}

	if(B == A) {
		fmpz_mat_t t;
		fmpz_mat_init(t, B->r, B->c);
		fmpz_mat_gram(t, A);
		fmpz_mat_swap_entrywise(B, t);
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
