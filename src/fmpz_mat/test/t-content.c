/*
    Copyright (C) 2015 Dharak Kharod

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_content, state)
{
	int i;
	fmpz_mat_t A,B;
	fmpz_t scalar, gcd_mat, temp;

	for (i = 0; i < 100 * flint_test_multiplier(); i++)
	{
		slong r, c;

		r = n_randint(state, 50);
		c = n_randint(state, 50);

		fmpz_mat_init(A, r, c);
		fmpz_mat_init(B, r, c);

		fmpz_init(scalar);
		fmpz_init(gcd_mat);
		fmpz_init(temp);

		fmpz_mat_randtest(B, state, 100);

		fmpz_mat_content(gcd_mat, B);

		if (r == 0 || c == 0)
		{
			if (fmpz_is_zero(gcd_mat))
			{
				goto cleanup;
			}
			else
			{
				flint_printf("FAIL!\n");
				flint_abort();
			}
		}

		fmpz_randtest_not_zero(scalar, state, 50);

		fmpz_mat_scalar_mul_fmpz(A, B, scalar);

		fmpz_mat_content(temp, A);

		fmpz_mul(gcd_mat, gcd_mat, scalar);

		if (fmpz_cmpabs(gcd_mat, temp) != 0)
		{
			flint_printf("FAIL!\n");
			flint_abort();
		}

cleanup:

		fmpz_mat_clear(A);
		fmpz_mat_clear(B);

		fmpz_clear(scalar);
		fmpz_clear(temp);
		fmpz_clear(gcd_mat);
	}

    TEST_FUNCTION_END(state);
}
