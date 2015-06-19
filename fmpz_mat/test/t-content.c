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

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "nmod_mat.h"
#include "ulong_extras.h"
#include "long_extras.h"


int main()
{
	int i;
	fmpz_mat_t A,B;
	fmpz_t scalar, gcd_mat, temp;
	FLINT_TEST_INIT(state);

	flint_printf("fmpz_mat_content....");
	fflush(stdout);

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
				continue;
			}
			else
			{
				flint_printf("FAIL!\n");
				abort();	
			}
		}
		
		fmpz_randtest_not_zero(scalar, state, 50);

		fmpz_mat_scalar_mul_fmpz(A, B, scalar);

		fmpz_mat_content(temp, A);
		
		fmpz_mul(gcd_mat, gcd_mat, scalar);

		if (fmpz_cmpabs(gcd_mat, temp) != 0)
		{
			flint_printf("FAIL!\n");
			abort();
		}

		fmpz_mat_clear(A);
		fmpz_mat_clear(B);
	
		fmpz_clear(scalar);
		fmpz_clear(temp);
		fmpz_clear(gcd_mat);
	}


	FLINT_TEST_CLEANUP(state);
	flint_printf("PASS\n");
	return 0;
}

