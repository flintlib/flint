/*
    Copyright (C) 2015 Dharak Kharod

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_sparse_mat.h"
#include "ulong_extras.h"
#include "long_extras.h"


int main()
{
	slong rep, r, c, min_nnz, max_nnz, bits, nreps = 100;
	fmpz_sparse_mat_t A,B;
	fmpz_t scalar, gcd_mat, temp;
	FLINT_TEST_INIT(state);

	flint_printf("fmpz_sparse_mat_content....");
	fflush(stdout);

	for (rep = 0; rep < nreps; ++rep)
	{
		r = n_randint(state, 50);
		c = n_randint(state, 50);
		min_nnz = 0;
		max_nnz = c;
		do bits = n_randint(state, 256);
		while (bits <= UWORD(1));

		fmpz_sparse_mat_init(A, r, c);
		fmpz_sparse_mat_init(B, r, c);

		fmpz_init(scalar);
		fmpz_init(gcd_mat);
		fmpz_init(temp);
		
		fmpz_sparse_mat_randtest(B, state, min_nnz, max_nnz, bits);
		
		fmpz_sparse_mat_content(gcd_mat, B);

		if (r == 0 || c == 0)
		{
			if (!fmpz_is_zero(gcd_mat))
			{
				flint_printf("FAIL!\n");
				abort();	
			}
		} else {	
			fmpz_randtest_not_zero(scalar, state, 50);

			fmpz_sparse_mat_scalar_mul_fmpz(A, B, scalar);

			fmpz_sparse_mat_content(temp, A);
			
			fmpz_mul(gcd_mat, gcd_mat, scalar);

			if (fmpz_cmpabs(gcd_mat, temp) != 0)
			{
				flint_printf("FAIL!\n");
				abort();
			}
		}
		fmpz_sparse_mat_clear(A);
		fmpz_sparse_mat_clear(B);
	
		fmpz_clear(scalar);
		fmpz_clear(temp);
		fmpz_clear(gcd_mat);
	}


	FLINT_TEST_CLEANUP(state);
	flint_printf("PASS\n");
	return 0;
}

