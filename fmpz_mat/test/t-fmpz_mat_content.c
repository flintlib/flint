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

    Copyright (C) 2015 Prabhdeep Singh Walia

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"
#include "long_extras.h"
int
main(void)
{
    int i;
    
    FLINT_TEST_INIT(state);

    flint_printf("fmpz_mat_content....");
    fflush(stdout);
    
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
	fmpz_mat_t X, Y;
	fmpz_t k, matgcd, tmp;
	
	slong rows, columns;
	rows = n_randint(state, 50);
	columns = n_randint(state, 50);

	fmpz_mat_init(Y, rows, columns);
	fmpz_mat_init(X, rows, columns);

	fmpz_init(k);
	fmpz_init(matgcd);
	fmpz_init(tmp);
		

	fmpz_mat_randtest(Y, state, 50);
	
	fmpz_mat_content(matgcd, Y);	
	fmpz_randtest_not_zero(k, state, 50);
	fmpz_mat_scalar_mul_fmpz(X, Y, k);
	fmpz_mat_content(tmp, X);
	fmpz_mul(matgcd, matgcd, k);

	if (fmpz_cmpabs(matgcd, tmp) != 0)
	{
		flint_printf("FAIL!\n");
		abort();
	}

	fmpz_mat_clear(X);
	fmpz_mat_clear(Y);

	fmpz_clear(k);
	fmpz_clear(tmp);
	fmpz_clear(matgcd);
      
    }


    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
