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

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpq.h"
#include "fmpq_mat.h"
#include "fmpq_vec.h"

int main(void)
{
	int i, result;
    FLINT_TEST_INIT(state);
    

    flint_printf("lll....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A, B, C;
        fmpq_mat_t Bq, Cq;
        fmpq_t alpha, onefourth, one;
        fmpz_t temp;

        slong m, n, bits;

		do {
			m = n_randint(state, 10);
			n = n_randint(state, 10);
		} while (m > n);

        bits = 1 + n_randint(state, 100);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, m, n);
        fmpz_mat_init(C, m, n);
        fmpq_mat_init(Bq, m, n);
        fmpq_mat_init(Cq, m, n);
        fmpq_init(alpha);
        fmpq_init(onefourth);
        fmpq_init(one);
        fmpz_init(temp);
		
		fmpz_mat_randrank(A, state, m, bits);
		fmpq_set_si(onefourth, 1, 4);
		fmpq_one(one);
		do {
			fmpq_randtest(alpha, state, bits);
		} while (fmpq_cmp(alpha, onefourth) <= 0 || fmpq_cmp(alpha, one) >= 0);
        
        fmpz_mat_rref(B, temp, A);
        fmpq_mat_set_fmpz_mat_div_fmpz(Bq, B, temp);
        
        fmpz_mat_lll(A, A, alpha);
        
		fmpz_mat_rref(C, temp, A);
		fmpq_mat_set_fmpz_mat_div_fmpz(Cq, C, temp);
        
        result = fmpq_mat_equal(Bq, Cq);
        
        if(!result) {
			flint_printf("FAIL:\n");
			flint_printf("B:\n");
			fmpz_mat_print_pretty(B);
			flint_printf("C:\n");
			fmpz_mat_print_pretty(C);
			abort();
		}
		
        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpq_mat_clear(Bq);
        fmpq_mat_clear(Cq);
        fmpq_clear(alpha);
        fmpq_clear(onefourth);
        fmpq_clear(one);
        fmpz_clear(temp);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
